# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# General pipeline parameters:

# Parent directory to pipeline workspace:
WORKSPACE=.
# Pipeline name:
PIPELINE_NAME=assembly-polish
# Pipeline working directory:
WDIR=$(WORKSPACE)/$(PIPELINE_NAME)
# Results directory:
RES=$(WDIR)/results
# Pipeline git repo:
REPO=https://github.com/nanoporetech/ont-assembly-polish.git

# Custom pipeline parameters:

# Input files:

 
NANOPORE_READS=
ILLUMINA_READS_PAIR1=
ILLUMINA_READS_PAIR2=

# Canu configuration:
CANU_GENOME_SIZE=12m
# Extra options passed to canu, refer to http://canu.readthedocs.io for more information:
CANU_PARAMETERS=

# Racon configuration:
# Use racon for polishing or not:
USE_RACON=yes
# Extra options passed to racon:
RACON_PARAMETERS=

# BWA configuration:
# Extra options passed to BWA mem:
BWA_PARAMETERS=

# Pilon configuration:
# Use pilon for polishing or not:
USE_PILON=yes
# Extra parameters passed to pilon, refer to https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage
PILON_PARAMETERS=
# Location of pilon release to use:
PILON_URL=https://github.com/broadinstitute/pilon/releases/download/v1.20/pilon-1.20.jar
# Name of pilon jar file:
PILON_JAR=pilon-1.20.jar
# Maximum amount of memory allocated by pilon. Increase this if pilon fails because of memory issues.
PILON_MAX_MEM=32G

# Number of cores to use for multithreaded applications:
CORES=32
