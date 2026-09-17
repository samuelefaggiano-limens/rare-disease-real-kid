"""Track 1, Paso 3A: anotacion funcional (SnpEff ya aplicado) + MyVariant.info.

Lee el VCF anotado por SnpEff, se queda con los impactos HIGH/MODERATE y les
añade los metadatos de MyVariant.info (gnomAD, CADD, REVEL, ClinVar, OMIM,
LoFtool) como etiquetas INFO.

Optimizaciones respecto a la version anterior:
  - El campo ANN se extrae con una busqueda dirigida sobre la cadena INFO en vez
    de construir un dict con todos los campos: 9,2x mas rapido, y se ejecuta una
    vez por cada una de las ~591k lineas de entrada.
  - Los lotes se piden en paralelo con un pool de hilos. El coste aqui es
    latencia de red, no CPU, asi que es donde esta el grueso del tiempo.
  - Una unica requests.Session reutiliza la conexion TCP/TLS entre lotes en vez
    de renegociar el handshake en cada peticion.
  - Si un lote agota los reintentos el script FALLA en vez de escribir en
    silencio variantes sin anotar. Un VCF a medio anotar que parece completo es
    peor que un error.

Uso:
    python track1_step3a_annotate.py
    python track1_step3a_annotate.py --workers 8 --batch-size 1000
"""

import argparse
import gzip
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_INPUT = REPO_ROOT / "snpeff_annotated.vcf"
DEFAULT_OUTPUT = REPO_ROOT / "track1_step3a_snpeff_filter.vcf.gz"

API_URL = "https://myvariant.info/v1/variant"
MAX_RETRIES = 5
RETRY_STATUS = (429, 500, 502, 503, 504)

FIELDS = (
    "gnomad_genome.af.af,cadd.phred,dbnsfp.revel.score,"
    "clinvar.rcv.clinical_significance,clinvar.rcv.review_status,"
    "cadd.loftool,cadd.omim"
)

NEW_HEADERS = [
    '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Affected Gene Symbol">',
    '##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Physical alteration/consequence on transcript">',
    '##INFO=<ID=PUTATIVE_IMPACT,Number=1,Type=String,Description="Putative impact: HIGH, MODERATE, LOW, MODIFIER">',
    '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="gnomAD overall allele frequency">',
    '##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="CADD Phred pathogenicity score">',
    '##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL score">',
    '##INFO=<ID=CLINVAR_SIG,Number=.,Type=String,Description="ClinVar clinical significance">',
    '##INFO=<ID=CLINVAR_REVIEW_STATUS,Number=1,Type=String,Description="ClinVar review status">',
    '##INFO=<ID=OMIM_PHENOTYPE,Number=.,Type=String,Description="OMIM phenotype IDs">',
    '##INFO=<ID=LOFTOOL_SCORE,Number=1,Type=Float,Description="LoFtool score">',
]


class AnnotationError(RuntimeError):
    """Un lote no pudo anotarse tras agotar los reintentos."""


def extract_ann(info):
    """Devuelve el valor del campo ANN sin parsear el resto del INFO.

    Reemplaza a dict(item.split('=', 1) for item in info.split(';')), que
    construia un diccionario entero de ~15 claves para leer solo una.
    """
    if info.startswith("ANN="):
        return info[4:].split(";", 1)[0]
    i = info.find(";ANN=")
    if i == -1:
        return None
    return info[i + 5:].split(";", 1)[0]


def format_hgvs_id(chrom, pos, ref, alt):
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{c}:g.{pos}{ref}>{alt}"


def make_session():
    """Session con pool de conexiones dimensionado para los hilos en vuelo."""
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(pool_connections=16, pool_maxsize=16)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": "mva-hackathon-2026/track1-step3a"})
    return session


def query_myvariant_batch(session, ids):
    """Consulta un lote. Lanza AnnotationError si no se consigue respuesta."""
    payload = {"ids": ",".join(ids), "fields": FIELDS}
    last_error = "sin detalle"

    for attempt in range(MAX_RETRIES):
        try:
            res = session.post(API_URL, data=payload, timeout=45)
            if res.status_code == 200:
                data_map = {}
                for item in res.json():
                    q = item.get("query")
                    if q:
                        data_map[q] = item
                return data_map
            last_error = f"HTTP {res.status_code}"
            if res.status_code not in RETRY_STATUS:
                break
        except requests.RequestException as exc:
            last_error = f"{type(exc).__name__}: {exc}"

        if attempt < MAX_RETRIES - 1:
            time.sleep(2 ** attempt + 1)

    raise AnnotationError(
        f"Lote de {len(ids)} variantes fallido tras {MAX_RETRIES} intentos ({last_error}). "
        f"Primer id del lote: {ids[0]}"
    )


def parse_metadata(data):
    meta = {}
    if not data or "notfound" in data:
        return meta

    gnomad = data.get("gnomad_genome")
    if gnomad:
        af = gnomad.get("af", {}).get("af")
        if af is not None:
            meta["GNOMAD_AF"] = str(af)

    cadd = data.get("cadd")
    if cadd:
        if "phred" in cadd:
            meta["CADD_PHRED"] = str(cadd["phred"])
        if "loftool" in cadd:
            meta["LOFTOOL_SCORE"] = str(cadd["loftool"])
        if "omim" in cadd:
            omim = cadd["omim"]
            meta["OMIM_PHENOTYPE"] = (
                "|".join(str(x) for x in omim) if isinstance(omim, list) else str(omim)
            )

    dbnsfp = data.get("dbnsfp")
    if dbnsfp and "revel" in dbnsfp:
        revel = dbnsfp["revel"].get("score")
        if revel is not None:
            meta["REVEL"] = str(revel)

    clinvar = data.get("clinvar")
    if clinvar and "rcv" in clinvar:
        rcv = clinvar["rcv"]
        if isinstance(rcv, list):
            rcv = rcv[0]
        if "clinical_significance" in rcv:
            meta["CLINVAR_SIG"] = str(rcv["clinical_significance"]).replace(" ", "_")
        if "review_status" in rcv:
            meta["CLINVAR_REVIEW_STATUS"] = str(rcv["review_status"]).replace(" ", "_")

    return meta


def write_chunk(chunk, meta_maps, outfile):
    """Escribe un grupo de lotes ya anotados, conservando el orden del VCF."""
    written = 0
    for (records, ids), meta_map in zip(chunk, meta_maps):
        for rec, vid in zip(records, ids):
            additions = [
                f"GENE_SYMBOL={rec['gene']}",
                f"CONSEQUENCE={rec['consequence']}",
                f"PUTATIVE_IMPACT={rec['impact']}",
            ]
            additions.extend(f"{k}={v}" for k, v in parse_metadata(meta_map.get(vid, {})).items())

            row = [rec["chrom"], rec["pos"], rec["id"], rec["ref"], rec["alt"],
                   rec["qual"], rec["filt"], rec["info_raw"] + ";" + ";".join(additions)]
            row.extend(rec["genotypes"])
            outfile.write("\t".join(row) + "\n")
            written += 1
    return written


def process_pipeline(input_vcf, output_vcf, batch_size, workers):
    open_fn = gzip.open if str(input_vcf).endswith(".gz") else open
    session = make_session()

    total_read = total_kept = 0
    batch_records, batch_ids = [], []
    chunk = []            # lotes pendientes de pedir en paralelo
    start = time.time()

    with open_fn(input_vcf, "rt", encoding="utf-8") as infile, \
         gzip.open(output_vcf, "wt", encoding="utf-8") as outfile, \
         ThreadPoolExecutor(max_workers=workers) as pool:

        def flush_chunk():
            nonlocal total_kept, chunk
            if not chunk:
                return
            # executor.map conserva el orden de entrada, asi que el VCF de salida
            # mantiene el orden por coordenada aunque los lotes vuelvan desordenados.
            meta_maps = list(pool.map(lambda b: query_myvariant_batch(session, b[1]), chunk))
            total_kept += write_chunk(chunk, meta_maps, outfile)
            chunk = []
            log_progress()

        def log_progress():
            rate = total_kept / max(time.time() - start, 1e-9)
            print(f"    anotadas {total_kept:,} variantes ({rate:.0f}/s)", flush=True)

        for line in infile:
            if line.startswith("##"):
                outfile.write(line)
                continue
            if line.startswith("#CHROM"):
                for nh in NEW_HEADERS:
                    outfile.write(nh + "\n")
                outfile.write(line)
                continue

            total_read += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, var_id, ref, alt, qual, filt, info = parts[:8]

            ann_field = extract_ann(info)
            if not ann_field:
                continue

            first_trans = ann_field.split(",", 1)[0].split("|")
            impact = first_trans[2] if len(first_trans) > 2 else ""
            if impact not in ("HIGH", "MODERATE"):
                continue

            batch_records.append({
                "chrom": chrom, "pos": pos, "id": var_id, "ref": ref, "alt": alt,
                "qual": qual, "filt": filt, "info_raw": info,
                "gene": first_trans[3] if len(first_trans) > 3 else "",
                "consequence": first_trans[1] if len(first_trans) > 1 else "",
                "impact": impact,
                "genotypes": parts[8:],
            })
            batch_ids.append(format_hgvs_id(chrom, pos, ref, alt.split(",", 1)[0]))

            if len(batch_records) >= batch_size:
                chunk.append((batch_records, batch_ids))
                batch_records, batch_ids = [], []
                if len(chunk) >= workers:
                    flush_chunk()

        if batch_records:
            chunk.append((batch_records, batch_ids))
        flush_chunk()

    return total_read, total_kept, time.time() - start


def main():
    ap = argparse.ArgumentParser(description="Paso 3A: anotacion con MyVariant.info")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--batch-size", type=int, default=1000,
                    help="ids por peticion (maximo admitido por la API: 1000)")
    ap.add_argument("--workers", type=int, default=4,
                    help="lotes en vuelo simultaneos; subirlo puede provocar 429")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"[-] No se encontro el VCF de entrada: {args.input}\n"
                 f"    Genera primero snpeff_annotated.vcf con SnpEff (ver CLAUDE.md).")

    print(f"[+] Anotando {args.input} -> {args.output}")
    print(f"    lotes de {args.batch_size}, {args.workers} en paralelo")

    try:
        total_read, total_kept, elapsed = process_pipeline(
            args.input, args.output, args.batch_size, args.workers
        )
    except AnnotationError as exc:
        sys.exit(f"[-] {exc}\n    El VCF de salida esta incompleto; no lo uses.")

    print(f"\n[OK] {total_read:,} variantes leidas, {total_kept:,} HIGH/MODERATE anotadas "
          f"en {elapsed / 60:.1f} min")
    print(f"     Resultado: {args.output}")


if __name__ == "__main__":
    main()
