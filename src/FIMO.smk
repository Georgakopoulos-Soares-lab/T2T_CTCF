import os
from pathlib import Path
from collections import defaultdict

# Define input directories
motif_dir = "../data/motifs"
data_dir = "../data"

# Auto-detect all motif files
motifs, = glob_wildcards(f"{motif_dir}/{{motif}}.meme")

# Find all FASTA files recursively
FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}

# First pass: group fastas by parent directory
dir_to_fastas = defaultdict(list)
for path in Path(data_dir).rglob("*"):
    if path.suffix in FASTA_EXTENSIONS:
        dir_to_fastas[path.parent].append(path)

# Second pass: build keys, keeping stem only when multiple fastas share a dir
fasta_files = {}
for parent, paths in dir_to_fastas.items():
    for path in paths:
        rel_parent = parent.relative_to(data_dir)
        if len(paths) > 1:
            # Multiple fastas (mat/pat) -> keep stem to distinguish
            key = str(rel_parent / path.stem)
        else:
            # Single fasta -> just use the directory path
            key = str(rel_parent)
        fasta_files[key] = str(path)

fastas = list(fasta_files.keys())

def get_fasta(wildcards):
    return fasta_files[wildcards.fasta]

rule all:
    input:
        expand("../res/fimo_res/{fasta}/{motif}.tsv", fasta=fastas, motif=motifs),
        expand("../data/fimo_backround_models/{fasta}.bg", fasta=fastas)

# rule hg38_only:
#     input:
#         expand("../res/fimo_res/{fasta}/{motif}.tsv",
#                fasta=[f for f in fastas if "hg38" in f], motif=motifs),
#         expand("../data/fimo_backround_models/{fasta}.bg",
#                fasta=[f for f in fastas if "hg38" in f])

rule create_background:
    input:
        genome = get_fasta
    output:
        background_file = "../data/fimo_backround_models/{fasta}.bg"
    wildcard_constraints:
        fasta = ".+"
    params:
        outdir = lambda wildcards, output: str(Path(output.background_file).parent)
    shell:
        """
        mkdir -p {params.outdir}
        fasta-get-markov -m 1 {input.genome} {output.background_file}
        """

rule run_fimo:
    input:
        motif = f"{motif_dir}/{{motif}}.meme",
        genome = get_fasta,
        background_file = "../data/fimo_backround_models/{fasta}.bg",
        fimo_script = "FIMO.sh"
    output:
        tsv = "../res/fimo_res/{fasta}/{motif}.tsv"
    wildcard_constraints:
        fasta = ".+",
        motif = "[^/]+"
    params:
        outdir = lambda wildcards, output: str(Path(output.tsv).parent)
    shell:
        """
        mkdir -p {params.outdir}
        bash {input.fimo_script} {input.motif} {input.background_file} {input.genome} {params.outdir}
        """