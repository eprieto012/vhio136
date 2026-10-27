# Installation required

cutadapt only works in a venv: https://cutadapt.readthedocs.io/en/stable/installation.html

1. generate my_virtual_env in the working directory:

cd <path to the folder>
python3 -m venv my_virtual_env
source my_virtual_env/bin/activate
python3 -m venv ~/cutadapt-venv
~/cutadapt-venv/bin/python3 -m pip install --upgrade pip
~/cutadapt-venv/bin/pip install git+https://github.com/marcelm/cutadapt.git

2. test the installation:

python3 -m cutadapt --version
which cutadapt

# TCR library screen analysis pipeline

3. run main.py

cd <path to the folder>
source my_virtual_env/bin/activate
parallel --link -a samples1.r1.txt -a samples1.r2.txt python3 main.py data/{1} data/{2} TCR_lib_RIO136_AAV.fasta -t 8

# Steps included in main.py

- Step 1: Fastp: trimming + qc + UMI extraction + dedup
- Step 2: umi-tools 
- Step 3: Fastp: trimming
- Step 4: Bowtie2 read alignment to TCR library