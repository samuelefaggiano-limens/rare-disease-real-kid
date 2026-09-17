"""Track 1, Paso 2: prefiltrado del VCF del paciente.

QC intrinseco (FILTER/GQ/DP/GT/AB) + sustraccion de variantes comunes contra
dbSNP Common. Mismos umbrales que track1_step2_exp01-vcf_prefiltering.ipynb.

Estrategia: tres pasadas SECUENCIALES, sin busquedas aleatorias y con memoria
acotada (~250 MB), en lugar de una consulta tabix por variante.

  Pasada 1  VCF paciente -> VCF temporal con los que superan el QC intrinseco.
            Guarda en memoria la clave (chrom,pos,ref,alt) de cada superviviente.
  Pasada 2  Recorre dbSNP una vez en binario y descarta del conjunto toda clave
            que aparezca alli. Lo que queda es el conjunto de variantes raras.
  Pasada 3  Reescribe el VCF temporal quedandose solo con las raras.

Por que asi y no de otra forma:
  - Cargar dbSNP entero en un set de tuplas cuesta ~287 bytes por variante,
    unos 11 GB para sus ~38M de registros: inviable.
  - Una consulta tabix por variante implica construir ~2,6M de iteradores y
    parsear cada registro de dbSNP como objeto completo. Aqui dbSNP se lee como
    bytes crudos partiendo solo los 5 primeros campos, sin tocar su INFO.
  - Las claves se guardan exactas, no hasheadas: una colision podria borrar la
    variante causal, y ese riesgo no compensa ahorrar ~100 MB.

Uso:
    python run_step2_prefilter.py
    python run_step2_prefilter.py --input ruta.vcf.gz --dbsnp common.vcf.gz -o out.vcf.gz
"""

import argparse
import gzip
import os
import sys
import time
from pathlib import Path

import cyvcf2

REPO_ROOT = Path(__file__).resolve().parents[2]

MIN_GQ = 20
MIN_DP = 10
HET_AB_MIN, HET_AB_MAX = 0.30, 0.70
HOM_ALT_AB_MIN = 0.90


def log(msg):
    print(msg, flush=True)


def variant_keys(chrom, pos, ref, alts):
    """Claves exactas de una variante, una por alelo alternativo.

    Normaliza el prefijo 'chr' porque el VCF del paciente y dbSNP usan contigs
    numericos, pero conviene no depender de ello.
    """
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return [f"{c}:{pos}:{ref}:{a}" for a in alts]


def pass1_intrinsic_qc(input_vcf, temp_vcf):
    """QC intrinseco. Devuelve el conjunto de claves supervivientes."""
    vcf_in = cyvcf2.VCF(input_vcf)
    vcf_in.add_info_to_header({
        "ID": "AB_FILTER",
        "Description": "Allelic Balance for Heterozygous calls",
        "Type": "Float",
        "Number": "1",
    })
    vcf_out = cyvcf2.Writer(temp_vcf, vcf_in)

    # Contadores locales: 3,4x mas rapidos que stats['clave'] += 1 en un bucle
    # de 5M iteraciones.
    total = dropped_filter = dropped_missing = 0
    dropped_gq = dropped_dp = dropped_gt = dropped_ab = 0
    survivors = set()
    start = time.time()

    for variant in vcf_in:
        total += 1

        # Los filtros se aplican de mas barato a mas caro. variant.format() crea
        # un array numpy en cada llamada, asi que AD no se pide hasta que GQ, DP
        # y GT ya han descartado lo que tenian que descartar.
        if variant.FILTER is not None and variant.FILTER != "PASS":
            dropped_filter += 1
            continue

        gq_array = variant.format("GQ")
        if gq_array is None:
            dropped_missing += 1
            continue
        if gq_array[0][0] < MIN_GQ:
            dropped_gq += 1
            continue

        dp_array = variant.format("DP")
        if dp_array is None:
            dropped_missing += 1
            continue
        if dp_array[0][0] < MIN_DP:
            dropped_dp += 1
            continue

        gt_type = variant.gt_types[0]
        if gt_type == 0 or gt_type == 2:  # HOM_REF y UNKNOWN
            dropped_gt += 1
            continue

        ad_array = variant.format("AD")
        if ad_array is None:
            dropped_missing += 1
            continue

        ref_depth = ad_array[0][0]
        alt_depth = ad_array[0][1] if len(ad_array[0]) > 1 else 0
        calc_dp = ref_depth + alt_depth
        if calc_dp == 0:
            dropped_ab += 1
            continue

        ab = alt_depth / calc_dp
        if gt_type == 1:
            if ab < HET_AB_MIN or ab > HET_AB_MAX:
                dropped_ab += 1
                continue
        elif gt_type == 3:
            if ab < HOM_ALT_AB_MIN:
                dropped_ab += 1
                continue

        variant.INFO["AB_FILTER"] = round(float(ab), 3)
        vcf_out.write_record(variant)
        survivors.update(variant_keys(variant.CHROM, variant.POS, variant.REF, variant.ALT))

        if total % 500_000 == 0:
            log(f"    {total:,} variantes leidas | {len(survivors):,} supervivientes "
                f"| {time.time() - start:.0f}s")

    vcf_in.close()
    vcf_out.close()

    stats = {
        "total": total,
        "dropped_filter": dropped_filter,
        "dropped_missing_metrics": dropped_missing,
        "dropped_gq": dropped_gq,
        "dropped_dp": dropped_dp,
        "dropped_gt": dropped_gt,
        "dropped_ab": dropped_ab,
    }
    return survivors, stats


def pass2_subtract_common(dbsnp_vcf, survivors):
    """Recorre dbSNP una vez y elimina del conjunto lo que sea comun.

    Lee en binario y parte solo los 5 primeros campos: el INFO de dbSNP es
    enorme y no se necesita para nada.
    """
    before = len(survivors)
    discard = survivors.discard  # se resuelve una vez, no 38M de veces
    lines = 0
    start = time.time()

    with gzip.open(dbsnp_vcf, "rb") as f:
        for line in f:
            lines += 1
            if line[0] == 35:  # '#'
                continue

            parts = line.split(b"\t", 5)
            if len(parts) < 5:
                continue

            chrom = parts[0].decode("ascii", "ignore")
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            pos = parts[1].decode("ascii", "ignore")
            ref = parts[3].decode("ascii", "ignore")
            alts = parts[4].decode("ascii", "ignore")

            prefix = f"{chrom}:{pos}:{ref}:"
            if "," in alts:
                for alt in alts.split(","):
                    discard(prefix + alt)
            else:
                discard(prefix + alts)

            if lines % 5_000_000 == 0:
                log(f"    {lines:,} lineas de dbSNP | {len(survivors):,} siguen siendo raras "
                    f"| {time.time() - start:.0f}s")

    return before - len(survivors)


def pass3_write_rare(temp_vcf, output_vcf, survivors):
    """Reescribe el temporal quedandose solo con las variantes raras."""
    vcf_in = cyvcf2.VCF(temp_vcf)
    vcf_out = cyvcf2.Writer(output_vcf, vcf_in)
    kept = 0

    for variant in vcf_in:
        keys = variant_keys(variant.CHROM, variant.POS, variant.REF, variant.ALT)
        # Todas sus claves entraron en el conjunto en la pasada 1. Si falta
        # alguna es que la pasada 2 la elimino por comun.
        if any(k not in survivors for k in keys):
            continue
        vcf_out.write_record(variant)
        kept += 1

    vcf_in.close()
    vcf_out.close()
    return kept


def main():
    ap = argparse.ArgumentParser(description="Paso 2: prefiltrado del VCF (QC + variantes raras)")
    ap.add_argument("--input", default=str(REPO_ROOT / "data" / "WGS_EX2312012_HGWCNDSX7.vcf.gz"))
    ap.add_argument("--dbsnp", default=str(REPO_ROOT / "data" / "common_all_20180418.vcf.gz"))
    ap.add_argument("-o", "--output", default=str(REPO_ROOT / "data" / "track1_step2_rare_variants.vcf.gz"))
    ap.add_argument("--keep-temp", action="store_true", help="no borrar el VCF intermedio")
    args = ap.parse_args()

    for path, label in ((args.input, "VCF del paciente"), (args.dbsnp, "dbSNP Common")):
        if not os.path.exists(path):
            sys.exit(f"[-] No se encontro el {label}: {path}")

    temp_vcf = args.output.replace(".vcf.gz", "") + ".intrinsic_qc.tmp.vcf.gz"
    t0 = time.time()

    log(f"[1/3] QC intrinseco sobre {args.input}")
    survivors, stats = pass1_intrinsic_qc(args.input, temp_vcf)
    passed_intrinsic = stats["total"] - sum(v for k, v in stats.items() if k.startswith("dropped"))
    log(f"      {passed_intrinsic:,} de {stats['total']:,} superan el QC "
        f"({time.time() - t0:.0f}s)")

    log(f"[2/3] Sustraccion de comunes contra {args.dbsnp}")
    dropped_common = pass2_subtract_common(args.dbsnp, survivors)
    log(f"      {dropped_common:,} claves descartadas por comunes ({time.time() - t0:.0f}s)")

    log(f"[3/3] Escribiendo variantes raras en {args.output}")
    kept = pass3_write_rare(temp_vcf, args.output, survivors)

    if not args.keep_temp:
        for suffix in ("", ".tbi", ".csi"):
            try:
                os.remove(temp_vcf + suffix)
            except OSError:
                pass

    elapsed = time.time() - t0
    log("\n--- ESTADISTICAS DE FILTRADO ---")
    log(f"Total de variantes:            {stats['total']:,}")
    log(f"Descartadas por FILTER:        {stats['dropped_filter']:,}")
    log(f"Descartadas sin metricas:      {stats['dropped_missing_metrics']:,}")
    log(f"Descartadas por GQ < {MIN_GQ}:       {stats['dropped_gq']:,}")
    log(f"Descartadas por DP < {MIN_DP}:       {stats['dropped_dp']:,}")
    log(f"Descartadas por GT (ref/unk):  {stats['dropped_gt']:,}")
    log(f"Descartadas por balance alel.: {stats['dropped_ab']:,}")
    log("--------------------------------")
    log(f"Superan el QC intrinseco:      {passed_intrinsic:,}")
    log(f"Descartadas por comunes:       {dropped_common:,}")
    log(f"Variantes raras finales:       {kept:,} "
        f"({kept / stats['total'] * 100:.1f}% del total)")
    log(f"\nCompletado en {elapsed / 60:.1f} min -> {args.output}")


if __name__ == "__main__":
    main()
