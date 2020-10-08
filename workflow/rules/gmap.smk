rule gmap_build:
    """Create gmap database from sequence file"""
    output:
        db = __RESOURCES__ / "gmap/db/{genome}.db.ok"
    input:
        seq = get_reference
    conda:
        "../envs/gmap.yaml"
    cache: True
    log:
        "logs/gmap_build/{genome}.log"
    wrapper:
        "file:///home/peru/dev/snakemake-workflows/rna-seq-map-smk/workflow/wrappers/bio/gmap/build"


rule gmap_map:
    """Map transcriptome to genome database"""
    output:
        res = report("data/interim/gmap/{genome}-{sample}.psl", caption="../report/gmap.rst", category="Gmap mapping")
    input:
        db = "resources/gmap/db/{genome}.db.ok",
        sample = get_sample,
        log = "logs/gmap_build/{genome}.log"
    conda:
        "../envs/gmap.yaml"
    threads:
        1
    log:
        "logs/gmap_map/{genome}-{sample}.psl.log"
    wrapper:
        "file:///home/peru/dev/snakemake-workflows/rna-seq-map-smk/workflow/wrappers/bio/gmap/map"


rule gmap_unzip_fasta:
    """Input files to gmap must be unzipped"""
    input:
        zipped = "{prefix}{fa}.gz"
    output:
        unzipped = temp("{prefix}{fa}")
    log:
        "logs/gmap_unzip_fasta/{prefix}{fa}.log"
    shell:
        "gzip -vdc {input.zipped} 2> {log} > {output.unzipped}"
