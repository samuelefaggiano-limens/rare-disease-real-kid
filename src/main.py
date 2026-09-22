import os
import sys
import subprocess
import time
import argparse
from pathlib import Path

# Add src to sys.path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.download import ensure_data

def run_step(name, script_path, args=None):
    if args is None:
        args = []
    print(f"\n======================================")
    print(f"=== Running {name} ===")
    print(f"======================================")
    cmd = [sys.executable, str(script_path)] + args
    start = time.time()
    res = subprocess.run(cmd)
    if res.returncode != 0:
        print(f"[-] {name} failed with exit code {res.returncode}")
        sys.exit(res.returncode)
    print(f"[OK] {name} completed in {time.time() - start:.1f}s")

def run_snpeff(force=False):
    """Runs snpEff to generate the annotated VCF for Step 3A."""
    input_vcf = REPO_ROOT / "data" / "track1_step2_rare_variants.vcf.gz"
    output_vcf = REPO_ROOT / "data" / "snpeff_annotated.vcf"
    
    # If the AlphaGenome output or MyVariant output already exists, we might not even need this.
    final_step3_vcf = REPO_ROOT / "data" / "track1_step3_annotation_and_filtering.vcf.gz"
    step3a_vcf = REPO_ROOT / "data" / "track1_step3a_snpeff_filter.vcf.gz"
    
    if not force and (final_step3_vcf.exists() or step3a_vcf.exists() or output_vcf.exists()):
        print("\n[+] Downstream VCFs for Step 3 exist. Skipping snpEff run.")
        return
        
    print("\n======================================")
    print("=== Running snpEff ===")
    print("======================================")
    
    snpeff_jar = REPO_ROOT / "data" / "snpEff" / "snpEff.jar"
    cmd = ["java", "-Xmx8g", "-jar", str(snpeff_jar), "-canon", "-noStats", "GRCh38.99", str(input_vcf)]
    
    start = time.time()
    try:
        with open(output_vcf, "w") as f_out:
            res = subprocess.run(cmd, stdout=f_out)
        if res.returncode != 0:
            print(f"[-] snpEff failed with exit code {res.returncode}")
            os.remove(output_vcf)
            sys.exit(res.returncode)
        print(f"[OK] snpEff completed in {time.time() - start:.1f}s")
    except FileNotFoundError:
        print("[-] 'java' not found. Please install Java to run snpEff, or ensure intermediate files are present.")
        if output_vcf.exists():
            os.remove(output_vcf)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="End-to-end Track 1 pipeline.")
    parser.add_argument("--force", action="store_true", help="Force re-run of all steps")
    parser.add_argument("--username", default="model1", help="Username for final output")
    parser.add_argument("--model", default="track1_submission", help="Model name for final output")
    args = parser.parse_args()

    # Step 0: Ensure Data
    ensure_data()
    
    # Step 2: VCF Pre-filtering
    step2_out = REPO_ROOT / "data" / "track1_step2_rare_variants.vcf.gz"
    if args.force or not step2_out.exists():
        run_step("Step 2 (Pre-filtering)", REPO_ROOT / "src" / "step2_prefilter.py", [
            "--input", str(REPO_ROOT / "data" / "WGS_EX2312012_HGWCNDSX7.vcf.gz"),
            "--dbsnp", str(REPO_ROOT / "data" / "common_all_20180418.vcf.gz"),
            "--output", str(step2_out)
        ])
    else:
        print("\n[+] Step 2 output exists, skipping.")
    
    # Run snpEff before Step 3A
    run_snpeff(force=args.force)
    
    # Step 3A: Annotate with MyVariant.info
    step3a_out = REPO_ROOT / "data" / "track1_step3a_snpeff_filter.vcf.gz"
    step3_final = REPO_ROOT / "data" / "track1_step3_annotation_and_filtering.vcf.gz"
    if args.force or (not step3a_out.exists() and not step3_final.exists()):
        run_step("Step 3A (MyVariant.info Annotation)", REPO_ROOT / "src" / "step3a_annotate.py", [
            "--input", str(REPO_ROOT / "data" / "snpeff_annotated.vcf"),
            "--output", str(step3a_out)
        ])
    else:
        print("\n[+] Step 3A output exists, skipping.")
    
    # Step 3B: AlphaGenome
    if args.force or not step3_final.exists():
        run_step("Step 3B (AlphaGenome)", REPO_ROOT / "src" / "step3b_alphagenome.py", [
            "--input", str(step3a_out),
            "--output", str(step3_final)
        ])
    else:
        print("\n[+] Step 3B output exists, skipping.")
    
    # Step 4: Phenotype Matching & Prioritization
    step4_out = REPO_ROOT / "data" / "track1_step4_candidates.tsv"
    if args.force or not step4_out.exists():
        run_step("Step 4 (Phenotype Matching)", REPO_ROOT / "src" / "step4_prioritize.py", [
            "--input", str(step3_final),
            "--docx", str(REPO_ROOT / "data" / "Challenge_Clinical_Phenotype_1.docx"),
            "--obo", str(REPO_ROOT / "data" / "hp.obo"),
            "--g2p", str(REPO_ROOT / "data" / "genes_to_phenotype.txt"),
            "--candidates-out", str(step4_out),
            "--ranking-out", str(REPO_ROOT / "data" / "track1_step4_gene_ranking.csv")
        ])
    else:
        print("\n[+] Step 4 output exists, skipping.")
    
    # Step 5: Ranking and Submission
    run_step("Step 5 (Ranking and Submission)", REPO_ROOT / "src" / "step5_submission.py", [
        "--input", str(step4_out),
        "--username", args.username,
        "--model", args.model,
        "--max-rows", "10"
    ])
    
    print("\n[+] Pipeline completed successfully!")
    print(f"[+] Final submission file saved to: output/track1/{args.username}_{args.model}.csv")

if __name__ == "__main__":
    main()
