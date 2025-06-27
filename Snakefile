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

match config.get("environment", "m2"):
    case "scicore":
        CONSENSUS_CMD="cp ~/nextclade_gisaid/sars-cov-2/pre-processed/synthetic.fasta"
    case "m2":
        CONSENSUS_CMD="cp ~/code/nextclade_data_workflows/sars-cov-2/pre-processed/synthetic.fasta"
    case "scp":
        CONSENSUS_CMD="scp -q roemer0001@login-transfer.scicore.unibas.ch:~/nextclade_gisaid/sars-cov-2/pre-processed/synthetic.fasta"
    case _:
        raise ValueError(f"Unknown environment: {config.get('environment', 'scicore')}")

rule get_consensus:
    output:
        "build/pango-consensus-sequences_genome-nuc_unsorted.fasta",
    params:
        command=CONSENSUS_CMD,
    shell:
        # To be adjusted if repos move
        """
        {params.command} {output}
        """


rule sort_sequences:
    input:
        "build/pango-consensus-sequences_genome-nuc_unsorted.fasta",
    output:
        "build/pango-consensus-sequences_genome-nuc.fasta",
    shell:
        "seqkit sort -N {input} >{output}"


rule run_nextclade:
    input:
        "build/pango-consensus-sequences_genome-nuc.fasta",
    output:
        tsv="build/nextclade.tsv",
    shell:
        "nextclade run --in-order -d sars-cov-2 --output-all build {input}"


rule download_open_metadata:
    output:
        "build/open_metadata.tsv.zst",
    shell:
        """
        aws s3 cp --no-progress s3://nextstrain-data/files/ncov/open/metadata.tsv.zst {output}
        """


rule find_open_lineages:
    """
    Publish consensus sequences only if there are 3 or more open genomes
    """
    input:
        rules.download_open_metadata.output,
    output:
        "build/open_lineages.txt",
    shell:
        """
        zstdcat {input} | \
        tsv-summarize -H --group-by Nextclade_pango --count | \
        tsv-filter -H --ge 'count:3' | \
        tsv-select -H -f1 >{output}
        """


rule compress_nuc:
    input:
        fasta="build/pango-consensus-sequences_genome-nuc.fasta",
        open="build/open_lineages.txt",
    output:
        "data/pango-consensus-sequences_genome-nuc.fasta.zst",
    shell:
        """
        seqkit grep -f {input.open} {input.fasta} | \
        zstd -f --ultra -21 >{output}
        """


rule compress_translations:
    input:
        fasta="build/nextclade_gene_{gene}.translation.fasta",
        open="build/open_lineages.txt",
    output:
        "data/pango-consensus-sequences_{gene}-aa.fasta.zst",
    shell:
        """
        seqkit grep -f {input.open} {input.fasta} | \
        zstd -f --ultra -21 >{output}
        """


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
        "data/pango-consensus-sequences_genome-nuc.fasta.zst",
        expand("data/pango-consensus-sequences_{gene}-aa.fasta.zst", gene=genes),
    output:
        "commit.done",
    params:
        add="git add" if not config.get("test", False) else "#",
        commit="git commit -m 'Automatic data update'; git push"
        if not config.get("test", False)
        else "",
    shell:
        """
        {params.add} {input}
        git status
        {params.commit}
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
