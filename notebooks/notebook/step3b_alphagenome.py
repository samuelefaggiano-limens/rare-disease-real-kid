"""Track 1, Paso 3B: enriquecimiento regulatorio con AlphaGenome Atlas.

Añade AVI_SCORE, TARGET_GENE y REGULATORY_FEATURE a las variantes que salen del
Paso 3A. Es un paso de enriquecimiento puro: no descarta ninguna variante.

Optimizaciones respecto a la version anterior:
  - Una sola pasada sobre el VCF en vez de dos. La version anterior lo leia
    entero para construir la lista de registros y volvia a leerlo desde cero
    para escribir, teniendo ya todo en memoria.
  - Memoria acotada a un lote. Antes se guardaban en RAM los campos de todas las
    variantes a la vez.
  - La API se consulta por lotes en vez de mandar todas las variantes en una
    unica llamada, que con volumenes altos agota cuota o expira.
  - Si el cliente de Atlas no se puede inicializar, el script se detiene. Antes
    seguia adelante y escribia un VCF sin anotar que aparentaba estar completo.

Uso:
    python step3b_alphagenome.py
    python step3b_alphagenome.py --batch-size 500
"""

import argparse
import gzip
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_INPUT = REPO_ROOT / "track1_step3a_snpeff_filter.vcf.gz"
DEFAULT_OUTPUT = REPO_ROOT / "track1_step3b_regulatory_annotated.vcf.gz"

NEW_HEADERS = [
    '##INFO=<ID=AVI_SCORE,Number=1,Type=Float,Description="AlphaGenome Variant Impact pathogenicity and regulatory disruption score">',
    '##INFO=<ID=TARGET_GENE,Number=.,Type=String,Description="Target gene(s) predicted to be altered by the variant">',
    '##INFO=<ID=REGULATORY_FEATURE,Number=.,Type=String,Description="Regulatory feature or mechanism disrupted (e.g. Splicing, Enhancer, Promoter)">',
]


def format_variant_id(chrom, pos, ref, alt):
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{c}:{pos}:{ref}:{alt}"


def build_client():
    """Inicializa el cliente de Atlas o aborta con un mensaje util."""
    try:
        from alphagenome.atlas.atlas import AtlasClient
    except ImportError as exc:
        sys.exit(f"[-] No se pudo importar alphagenome ({exc}).\n"
                 f"    Instalalo con: pip install alphagenome")
    try:
        return AtlasClient()
    except Exception as exc:
        sys.exit(f"[-] No se pudo inicializar AtlasClient ({exc}).\n"
                 f"    Revisa la clave de API de Google DeepMind en el entorno.")


def query_atlas(client, variant_ids):
    """Consulta un lote y devuelve {id_variante: {campos}}."""
    results = client.query_variants(variant_ids, requested_scorers=["avi"])
    obs_df = getattr(results, "obs", None)
    if obs_df is None:
        return {}

    return {
        idx: {
            "AVI_SCORE": row.get("avi_score", ""),
            "TARGET_GENE": row.get("target_gene", ""),
            "REGULATORY_FEATURE": row.get("regulatory_feature", ""),
        }
        for idx, row in obs_df.iterrows()
    }


def flush_batch(client, batch, outfile):
    """Anota un lote y lo escribe. batch = [(parts, var_id), ...]."""
    atlas_data = query_atlas(client, [vid for _, vid in batch])

    for parts, var_id in batch:
        data = atlas_data.get(var_id)
        if data:
            additions = [f"{k}={v}" for k, v in data.items() if v != ""]
            if additions:
                parts[7] = parts[7] + ";" + ";".join(additions)
        outfile.write("\t".join(parts) + "\n")

    return len(batch)


def run_step3b(input_vcf, output_vcf, batch_size):
    client = build_client()
    print(f"[+] Enriqueciendo {input_vcf} -> {output_vcf}")

    batch = []
    total = 0
    header_done = False

    with gzip.open(input_vcf, "rt", encoding="utf-8") as infile, \
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
                total += flush_batch(client, batch, outfile)
                batch = []
                print(f"    {total:,} variantes enriquecidas", flush=True)

        if batch:
            total += flush_batch(client, batch, outfile)

    if not header_done:
        sys.exit(f"[-] {input_vcf} no tiene linea #CHROM; el VCF de entrada esta vacio o corrupto.")

    print(f"\n[OK] Paso 3B completado: {total:,} variantes enriquecidas -> {output_vcf}")
    return total


def main():
    ap = argparse.ArgumentParser(description="Paso 3B: enriquecimiento con AlphaGenome Atlas")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--batch-size", type=int, default=500,
                    help="variantes por llamada a la API")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"[-] No se encontro el VCF de entrada: {args.input}\n"
                 f"    Ejecuta antes el Paso 3A.")

    run_step3b(args.input, args.output, args.batch_size)


if __name__ == "__main__":
    main()
