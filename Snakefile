genes = [
    "S",
    "ORF1a",
    "ORF1b",
    "ORF3a",
    "E",
    "M",
    "N",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF9b",
]


rule all:
    input:
        commit_touchfile="commit.done",
        # "upload.done",
    shell:
        """
        rm {input.commit_touchfile}
        """


rule get_consensus:
    output:
        "build/pango-consensus-sequences_nuc_unsorted.fasta",
    shell:
        # To be adjusted if repos move
        "cp ~/nextclade_gisaid/sars-cov-2/pre-processed/synthetic.fasta {output}"


rule sort_sequences:
    input:
        "build/pango-consensus-sequences_nuc_unsorted.fasta",
    output:
        "build/pango-consensus-sequences_nuc.fasta",
    shell:
        "seqkit sort -N {input} >{output}"


rule run_nextclade:
    input:
        "build/pango-consensus-sequences_nuc.fasta",
    output:
        tsv="build/nextclade.tsv",
    shell:
        "nextclade run --in-order -d sars-cov-2 --output-all build {input}"


rule compress_nuc:
    input:
        "build/pango-consensus-sequences_{seq_type}.fasta",
    output:
        "data/pango-consensus-sequences_{seq_type}.fasta.zst",
    shell:
        "zstd -f --ultra -21 {input} >{output}"


rule compress_translations:
    input:
        "build/nextclade_gene_{gene}.translation.fasta",
    output:
        "data/pango-consensus-sequences_{gene}.fasta.zst",
    shell:
        "zstd -f --ultra -21 {input} >{output}"


rule create_json:
    input:
        nextclade_tsv="build/nextclade.tsv",
    output:
        "data/pango-consensus-sequences_summary.json",
    shell:
        """
        python scripts/create_json.py \
            --nextclade-tsv {input.nextclade_tsv} \
            --output {output}
        """


rule commit_results:
    input:
        "data/pango-consensus-sequences_summary.json",
        "data/pango-consensus-sequences_nuc.fasta.zst",
        expand("data/pango-consensus-sequences_{gene}.fasta.zst", gene=genes),
    output:
        "commit.done",
    shell:
        """
        git add {input}
        git status
        git commit -m "Update pango consensus sequences" || true
        git push
        touch {output}
        """


# rule upload_to_s3:
#     """
#     Upload of uncompressed files to S3
#     All fastas and json
#     """
#     input:
#     output:
#         "upload.done",
#     shell:
#         """
#         aws s3 sync data/ s3://nextstrain-data/files/pango-designation/
#         """
