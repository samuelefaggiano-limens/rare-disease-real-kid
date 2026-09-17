"""Track 1, Paso 3B: enriquecimiento regulatorio con AlphaGenome Atlas.

Añade a cada variante del Paso 3A las puntuaciones del scorer AVI (AlphaGenome
Variant Impact). Es un paso de enriquecimiento puro: no descarta ninguna
variante.

Etiquetas INFO que añade:
  AVI_SCORE           puntuacion cruda del modelo
  AVI_QUANTILE        cuantil calibrado (1 - CDF) de la variante frente a todos
                      los SNV del genoma
  AVI_PHRED           -10*log10(1 - cuantil). >= 20 es el 1 % superior del
                      genoma, >= 30 el 0,1 %
  REGULATORY_FEATURE  modalidad con mayor atribucion (MERGED_SPLICING,
                      ALPHAMISSENSE, MAX_ABS_RNA_SEQ, ...) del scorer
                      AVI_SCORE_FEATURE_IMPORTANCE
  AVI_STATUS          'ok' o el motivo por el que no se pudo puntuar
                      (p. ej. ref_mismatch cuando el alelo REF no coincide con
                      la referencia GRCh38 del servidor)

TARGET_GENE, previsto en docs/01-procedure.md, no lo devuelven los scorers AVI;
el gen sigue siendo el GENE_SYMBOL de snpEff del Paso 3A.

Como se usa la API (alphagenome 0.9):
  - El cliente se crea con atlas.create(api_key); AtlasClient() a secas no
    funciona (necesita el stub gRPC). La clave se lee de ALPHAGENOME_API_KEY
    en ~/.env (convencion de las skills de AlphaGenome) o del entorno.
  - query_variants() devuelve {scorer: AnnData}; obs['variant'] identifica la
    fila, X[:,0] es el score y layers['quantiles'] el cuantil.
  - El servidor valida el alelo REF: una sola variante mal referenciada hace
    fallar TODO el lote. Por eso, si un lote falla con ValueError, se reintenta
    variante a variante y solo las que fallan quedan marcadas en AVI_STATUS.
  - Los errores de infraestructura (red, cuota) se reintentan con backoff y,
    si persisten, el script ABORTA en vez de escribir un VCF a medio anotar.

Uso:
    python step3b_alphagenome.py
    python step3b_alphagenome.py --batch-size 100 --workers 8
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_INPUT = REPO_ROOT / "track1_step3a_snpeff_filter.vcf.gz"
DEFAULT_OUTPUT = REPO_ROOT / "track1_step3b_regulatory_annotated.vcf.gz"
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"]
MAX_RETRIES = 5

NEW_HEADERS = [
    '##INFO=<ID=AVI_SCORE,Number=1,Type=Float,Description="AlphaGenome Variant Impact raw score">',
    '##INFO=<ID=AVI_QUANTILE,Number=1,Type=Float,Description="AlphaGenome Variant Impact calibrated tail quantile (1-CDF) vs genome-wide SNVs">',
    '##INFO=<ID=AVI_PHRED,Number=1,Type=Float,Description="-10*log10(1-AVI_QUANTILE); >=20 is top 1 percent genome-wide">',
    '##INFO=<ID=REGULATORY_FEATURE,Number=1,Type=String,Description="Top contributing AlphaGenome modality (feature importance)">',
    '##INFO=<ID=AVI_STATUS,Number=1,Type=String,Description="ok, or reason the variant could not be scored">',
]


def log(msg):
    print(msg, flush=True)


def load_api_key() -> str:
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        env_path = os.path.expanduser("~/.env")
        if os.path.exists(env_path):
            with open(env_path) as fh:
                for line in fh:
                    if line.startswith("ALPHAGENOME_API_KEY="):
                        key = line.split("=", 1)[1].strip().strip('"').strip("'")
                        break
    if not key:
        sys.exit("[-] Falta ALPHAGENOME_API_KEY (en el entorno o en ~/.env).\n"
                 "    Guardala con: printf 'Enter ALPHAGENOME_API_KEY: ' && read -s v && echo \"ALPHAGENOME_API_KEY=$v\" >> ~/.env")
    return key


def build_client():
    try:
        from alphagenome.atlas import atlas
    except ImportError as exc:
        sys.exit(f"[-] No se pudo importar alphagenome ({exc}). Instalalo con: pip install alphagenome")
    try:
        return atlas.create(load_api_key())
    except Exception as exc:  # noqa: BLE001
        sys.exit(f"[-] No se pudo crear el cliente de Atlas ({exc}).")


def format_variant_id(chrom, pos, ref, alt):
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{c}:{pos}:{ref}:{alt}"


def to_atlas_str(vid: str) -> str:
    """chr15:100:C:T -> chr15:100:C>T (formato de genome.Variant.from_str)."""
    c, p, r, a = vid.split(":")
    return f"{c}:{p}:{r}>{a}"


def parse_results(res) -> dict[str, dict]:
    """{variant_str: {AVI_SCORE, AVI_QUANTILE, AVI_PHRED, REGULATORY_FEATURE}}."""
    out: dict[str, dict] = {}
    score_ad = res.get("AVI_SCORE")
    if score_ad is not None:
        X = score_ad.X.toarray() if hasattr(score_ad.X, "toarray") else score_ad.X
        q_layer = score_ad.layers.get("quantiles") if hasattr(score_ad, "layers") else None
        if q_layer is not None and hasattr(q_layer, "toarray"):
            q_layer = q_layer.toarray()
        for i, v in enumerate(score_ad.obs["variant"]):
            d = out.setdefault(str(v), {})
            d["AVI_SCORE"] = f"{float(X[i, 0]):.4f}"
            if q_layer is not None:
                q = float(q_layer[i, 0])
                d["AVI_QUANTILE"] = f"{q:.6f}"
                d["AVI_PHRED"] = f"{-10 * math.log10(max(1 - q, 1e-12)):.2f}"
    fi_ad = res.get("AVI_SCORE_FEATURE_IMPORTANCE")
    if fi_ad is not None:
        X = fi_ad.X.toarray() if hasattr(fi_ad.X, "toarray") else fi_ad.X
        names = list(fi_ad.var["name"])
        for i, v in enumerate(fi_ad.obs["variant"]):
            row = X[i]
            j = max(range(len(names)), key=lambda k: abs(float(row[k])))
            out.setdefault(str(v), {})["REGULATORY_FEATURE"] = names[j]
    return out


def query_batch(client, variant_ids: list[str], workers: int) -> tuple[dict[str, dict], dict[str, str]]:
    """Devuelve (anotaciones, errores_por_variante). Aborta si la API no responde."""
    from alphagenome.data import genome

    def run(ids):
        vs = [genome.Variant.from_str(to_atlas_str(v)) for v in ids]
        return client.query_variants(vs, requested_scorers=SCORERS, progress_bar=False, max_workers=workers)

    errors: dict[str, str] = {}
    for attempt in range(MAX_RETRIES):
        try:
            return parse_results(run(variant_ids)), errors
        except (ValueError, NotImplementedError) as exc:
            # Validacion del servidor (REF no coincide, contig desconocido,
            # tipo de variante no soportado): afecta a una variante, pero
            # tumba el lote. Se resuelve una a una.
            if len(variant_ids) == 1:
                errors[variant_ids[0]] = classify_error(exc)
                return {}, errors
            annotations: dict[str, dict] = {}
            for vid in variant_ids:
                a, e = query_batch(client, [vid], workers)
                annotations.update(a)
                errors.update(e)
            return annotations, errors
        except Exception as exc:  # noqa: BLE001  errores de red / gRPC / cuota
            wait = 2 ** attempt + 1
            log(f"    aviso: lote fallido ({type(exc).__name__}: {str(exc)[:120]}); reintento en {wait}s")
            time.sleep(wait)
    sys.exit(f"[-] La API de Atlas no respondio tras {MAX_RETRIES} intentos. "
             f"El VCF de salida esta incompleto; no lo uses.")


def classify_error(exc: Exception) -> str:
    msg = str(exc)
    if isinstance(exc, NotImplementedError):
        return "unsupported_by_atlas"
    if "reference base" in msg:
        return "ref_mismatch"
    if "chromosome" in msg.lower() or "contig" in msg.lower():
        return "unknown_contig"
    return "invalid:" + msg[:60].replace(";", ",").replace(" ", "_")


def is_snv(vid: str) -> bool:
    _, _, ref, alt = vid.split(":")
    return len(ref) == 1 and len(alt) == 1


def flush_batch(client, batch, outfile, workers, stats):
    # Atlas (alphagenome 0.9, 2026-09) responde NotImplementedError a todo
    # indel: se marcan sin gastar una llamada. Solo se consultan los SNV.
    ids = [vid for _, vid in batch if is_snv(vid)]
    annotations, errors = query_batch(client, ids, workers) if ids else ({}, {})
    for parts, vid in batch:
        if not is_snv(vid):
            errors[vid] = "indel_unsupported"
        data = annotations.get(to_atlas_str(vid)) or annotations.get(vid)
        additions = []
        if data:
            additions.extend(f"{k}={v}" for k, v in data.items())
            additions.append("AVI_STATUS=ok")
            stats["scored"] += 1
        else:
            additions.append(f"AVI_STATUS={errors.get(vid, 'no_result')}")
            stats["unscored"] += 1
        parts[7] = parts[7] + ";" + ";".join(additions)
        outfile.write("\t".join(parts) + "\n")
    return len(batch)


def run_step3b(input_vcf, output_vcf, batch_size, workers):
    client = build_client()
    log(f"[+] Enriqueciendo {input_vcf} -> {output_vcf}")

    batch, total, header_done = [], 0, False
    stats = {"scored": 0, "unscored": 0}
    start = time.time()
    open_in = gzip.open if str(input_vcf).endswith(".gz") else open

    with open_in(input_vcf, "rt", encoding="utf-8") as infile, \
         gzip.open(output_vcf, "wt", encoding="utf-8") as outfile:
        for line in infile:
            if line.startswith("##"):
                outfile.write(line)
                continue
            if line.startswith("#CHROM"):
                for nh in NEW_HEADERS:
                    outfile.write(nh + "\n")
                outfile.write(line)
                header_done = True
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[:5]
            batch.append((parts, format_variant_id(chrom, pos, ref, alt.split(",", 1)[0])))
            if len(batch) >= batch_size:
                total += flush_batch(client, batch, outfile, workers, stats)
                batch = []
                rate = total / max(time.time() - start, 1e-9)
                log(f"    {total:,} variantes ({rate:.1f}/s) | puntuadas {stats['scored']:,} | sin puntuar {stats['unscored']:,}")
        if batch:
            total += flush_batch(client, batch, outfile, workers, stats)

    if not header_done:
        sys.exit(f"[-] {input_vcf} no tiene linea #CHROM; el VCF de entrada esta vacio o corrupto.")

    log(f"\n[OK] Paso 3B: {total:,} variantes escritas, {stats['scored']:,} con AVI, "
        f"{stats['unscored']:,} sin puntuar (ver AVI_STATUS) en {(time.time() - start) / 60:.1f} min -> {output_vcf}")
    return total


def main():
    ap = argparse.ArgumentParser(description="Paso 3B: enriquecimiento con AlphaGenome Atlas (AVI)")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--batch-size", type=int, default=100, help="variantes por llamada")
    ap.add_argument("--workers", type=int, default=8, help="peticiones gRPC en vuelo dentro de un lote")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"[-] No se encontro el VCF de entrada: {args.input}\n    Ejecuta antes el Paso 3A.")
    run_step3b(args.input, args.output, args.batch_size, args.workers)


if __name__ == "__main__":
    main()
