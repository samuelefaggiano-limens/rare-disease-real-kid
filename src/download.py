import os
import sys
from pathlib import Path
from huggingface_hub import hf_hub_download, login

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"

REQUIRED_FILES = [
    "WGS_EX2312012_HGWCNDSX7.vcf.gz",
    "WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi",
    "Challenge_Clinical_Phenotype_1.docx",
    "common_all_20180418.vcf.gz",
    "common_all_20180418.vcf.gz.tbi",
    "hp.obo",
    "genes_to_phenotype.txt",
    "snpEff.zip"
]

def ensure_data():
    DATA_DIR.mkdir(exist_ok=True)
    repo_id = "SageBio/mva-hackathon-2026-data"
    
    missing_files = [f for f in REQUIRED_FILES if not (DATA_DIR / f).exists()]
    
    if not missing_files:
        print("[+] All required datasets are already present in data/")
        return

    print(f"[*] Missing files detected: {missing_files}")
    print("[*] Attempting to download from Hugging Face Hub...")
    
    # Check if we have token in env
    token = os.environ.get("HF_TOKEN")
    if token:
        login(token=token)
    
    for filename in missing_files:
        print(f"    Downloading {filename}...")
        try:
            downloaded_path = hf_hub_download(
                repo_id=repo_id,
                repo_type="dataset",
                filename=filename,
                local_dir=DATA_DIR
            )
            print(f"    [OK] Downloaded {filename}")
        except Exception as e:
            print(f"    [-] Failed to download {filename}: {e}")
            print("        Make sure you have access to the gated dataset and HF_TOKEN is set.")
            sys.exit(1)

if __name__ == "__main__":
    ensure_data()
