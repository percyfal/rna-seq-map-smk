__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

seq = snakemake.input.seq
touchdb = snakemake.output.db
zipflag = "-g" if re.search("(.gz|.gzip)$", snakemake.wildcards.genome) else ""

dirname = os.path.dirname(touchdb)
db = f"{snakemake.wildcards.genome}.db"

shell(
    "gmap_build {extra} "
    "{zipflag} "
    "-d {db} "
    "--dir {dirname} "
    "{seq} "
    "{log}"
)
shell(
    "touch {touchdb}"
)
