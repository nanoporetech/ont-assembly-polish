# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Pipeline targets:

# Nanopore assembly by canu:

CANU_PREFIX=canu
CANU_DIR=$(WDIR)/canu-assembly
CANU_CONTIGS=$(CANU_DIR)/$(CANU_PREFIX).contigs.fasta

canu_assembly: $(CANU_CONTIGS) $(WDT)
$(CANU_CONTIGS): $(NANOPORE_READS)
	@echo Assembling nanopore reads using canu.
	@canu\
		 -p $(CANU_PREFIX) \
		 -d $(CANU_DIR) \
		 genomeSize=$(CANU_GENOME_SIZE) \
		 -nanopore-raw $(NANOPORE_READS) $(CANU_PARAMETERS)

# Map nanopore reads to canu contings and use racon to perform correction based on nanopore reads only:

MINIMAP_OVERLAPS=$(WDIR)/minimap_overlaps.paf

ifeq ($(USE_RACON),yes)
RACON_CONTIGS=$(WDIR)/racon.contigs.fasta
else
RACON_CONTIGS=$(CANU_CONTIGS)
endif

# Racon has to be rebuilt on the current machine or it might fail.

build_racon: $(WDIR)/racon
$(WDIR)/racon: $(WDT)
	(cd $(WDIR) && git clone https://github.com/isovic/racon.git racon.build && cd racon.build && make modules && make tools && make -j) &&\
	 mv $(WDIR)/racon.build/bin/racon $(WDIR)/ && rm -fr racon.build

racon_correct: $(RACON_CONTIGS)
$(RACON_CONTIGS): $(NANOPORE_READS) $(CANU_CONTIGS) $(WDIR)/racon
ifeq ($(USE_RACON),yes)
	@echo Mapping nanopore reads onto canu contings using minimap.
	@minimap $(CANU_CONTIGS) $(NANOPORE_READS) > $(MINIMAP_OVERLAPS)
	@echo Correcting contigs using racon.
	@$(WDIR)/racon -t $(CORES) $(NANOPORE_READS) $(MINIMAP_OVERLAPS) $(CANU_CONTIGS) $(RACON_CONTIGS)
else
	@echo Skipping racon polishing.
endif

# Index contigs, map Illumina reads to contigs by BWA, sorting and indexing using samtools:

BWA_BAM_PREFIX=$(WDIR)/bwa_aligned_reads
BWA_BAM=$(BWA_BAM_PREFIX).bam

bwa_align: $(BWA_BAM)
$(BWA_BAM): $(RACON_CONTIGS) $(ILLUMINA_READS_PAIR1) $(ILLUMINA_READS_PAIR2)
ifeq ($(USE_PILON),yes)
	@echo Indexing contigs at $(RACON_CONTIGS).
	@bwa index $(RACON_CONTIGS)
	@echo Aligning Illumina reads using BWA mem.
	@bwa mem -t $(CORES) $(BWA_PARAMETERS) $(RACON_CONTIGS)  $(ILLUMINA_READS_PAIR1) $(ILLUMINA_READS_PAIR2)\
		| samtools view -S -b -u - | samtools sort - $(BWA_BAM_PREFIX)
	@samtools index $(BWA_BAM)
else
	@echo Skipping BWA alignment.
endif

# Correct contigs using pilon based on the Illumina reads:

PILON_CONTIGS=$(WDIR)/pilon.contigs.fasta

# Get pilon:
get_pilon: $(WDIR)/$(PILON_JAR)
$(WDIR)/$(PILON_JAR):
ifeq ($(USE_PILON),yes)
	(cd $(WDIR) && wget $(PILON_URL))
else
	@echo Skipping pilon download.
endif

pilon_correct: $(PILON_CONTIGS)
$(PILON_CONTIGS): $(RACON_CONTIGS) $(BWA_BAM) $(WDIR)/$(PILON_JAR)
ifeq ($(USE_PILON),yes)
	@echo Correcting contigs using pilon.
	@java -Xmx$(PILON_MAX_MEM) -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit -jar $(WDIR)/$(PILON_JAR) --threads $(CORES) --genome $(RACON_CONTIGS)\
		--bam $(BWA_BAM) --outdir $(WDIR) --output pilon.contigs $(PILON_PARAMETERS)
else
	@echo Skipping pilon polishing.
endif

all: $(PILON_CONTIGS)
	@echo
	@echo Analysis finishes.
	@echo Contigs assembled by Canu are at: $(CANU_CONTIGS)
	@echo Contigs polished by racon are at: $(RACON_CONTIGS)
	@echo Pilon corrected contigs are at: $(PILON_CONTIGS)

