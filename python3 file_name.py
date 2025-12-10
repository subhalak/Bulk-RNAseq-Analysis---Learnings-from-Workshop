import subprocess
import time

# List of all SRA IDs
sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540", "SRR7179541"
]

# 1) Download each .sra file
for sra_id in sra_numbers:
    print(f"\n=== Downloading: {sra_id} ===")
    prefetch_cmd = f"prefetch {sra_id}"
    print("Command:", prefetch_cmd)

    start = time.time()
    subprocess.call(prefetch_cmd, shell=True)
    end = time.time()

    print(f"⏱ Download time for {sra_id}: {(end-start)/60:.2f} minutes")


# 2) Convert each .sra file to fastq.gz
for sra_id in sra_numbers:
    sra_file = f"{sra_id}.sra"

    print(f"\n=== Generating FASTQ for: {sra_id} ===")
    fastq_dump_cmd = (
        f"fastq-dump --outdir fastq --gzip --skip-technical "
        f"--readids --read-filter pass --dumpbase --split-3 --clip {sra_file}"
    )

    print("Command:", fastq_dump_cmd)

    start = time.time()
    subprocess.call(fastq_dump_cmd, shell=True)
    end = time.time()

    print(f"⏱ FASTQ conversion time for {sra_id}: {(end-start)/60:.2f} minutes")

