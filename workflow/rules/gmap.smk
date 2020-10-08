rule gmap_build_local:
    """Create gmap genome database from local genome file"""
    output:
        db = "resources/gmap/db/{genome}.db.ok"
    input:
        seq = get_reference
    conda:
        "../envs/gmap.yaml"
    cache: True
    params:
        # NB: check in wrapper
        gunzip = lambda wildcards: "-g" if wildcards.genome.endswith(".gz") else "",
        db = lambda wildcards: f"{wildcards.genome}.db"
    log:
        "logs/gmap_build_local/{genome}.log"
    shell:
        """gmap_build {params.gunzip} -d {params.db} --dir $(dirname {output.db}) {input.seq} 2> {log} && touch {output.db}"""


rule gmap_map:
    """Map transcriptome to genome database"""
    output:
        res = report("data/interim/gmap/{genomedb}-{sample}.{gmap}.psl", caption="../report/gmap.rst", category="Gmap mapping")
    input:
        db = "resources/gmap/db/{genomedb}.db.ok",
        sample = get_sample
    conda:
        "../envs/gmap.yaml"
    wildcard_constraints:
        # NB: should be checked in wrapper
        gmap = "(gmap|gmapl)"
    params:
        # NB: check in wrapper? For now only allow psl
        format = 1,
        dir = lambda wildcards, input: str(Path(input.db).parent),
        db = lambda wildcards, input: str(Path(input.db).name).rstrip(".ok")
    threads:
        1
    log:
        "logs/gmap_map/{genomedb}-{sample}.{gmap}.psl.log"
    shell:
        "{wildcards.gmap} -t {threads} --dir {params.dir} --db {params.db} -f {params.format} {input.sample} > {output.res} 2> {log}"


rule gmap_unzip_fasta:
    """Input files to gmap must be unzipped"""
    input:
        zipped = "{prefix}.gz"
    output:
        unzipped = temp("{prefix}")
    log:
        "logs/gmap_unzip_fasta/{prefix}.log"
    shell:
        "gzip -vdc {input.zipped} 2> {log} > {output.unzipped}"
