ONT assembly and Illumina polishing pipeline
=============================================

This pipeline performs the following steps:
- Assembly of nanopore reads using [Canu](http://canu.readthedocs.io).
- Polish canu contigs using [racon](https://github.com/isovic/racon) (*optional*).
- Map a paired-end Illumina dataset onto the contigs obtained in the previous steps using [BWA](http://bio-bwa.sourceforge.net) mem.
- Perform correction of contigs using [pilon](https://github.com/broadinstitute/pilon/wiki) and the Illumina dataset.

Usage
-----

Edit `config.mk` to set input files and parameters. Specifying the following is mandatory:
- `NANOPORE_READS` - input nanopore reads (note that this **must** be a single valid fastq file, see [here](https://www.biostars.org/p/81924/) how to combine fastq files).
- `ILLUMINA_READS_PAIR1` - fastq with the first reads of the paired-end Illumina dataset.
- `ILLUMINA_READS_PAIR2` - fastq with the second reads of the paired-end Illumina dataset.
- `CANU_GENOME_SIZE` - genome size parameter passed to canu.
- `PILON_MAX_MEM` - maximum amount of memory used by pilon. Increase this is if pilon crashes because of running out of memory.

The number of cores used can be specified by `CORES` (set this to the number of CPUs in your machine).
Racon corrections can be disabled by setting `USE_RACON=no`. Pilon polishing can be disabled by setting `USE_PILON=no`.

Then issue issue `make all` to run the pipeline. Issue `make help` for a list of utility make targets. Issue `make clear_wdir` to delete the working directory (including all results!).

Using through docker
--------------------

The easiest way to use the pipeline is through docker. First [install docker](https://docs.docker.com/engine/installation/), then issue the following to build the
container:

```bash
cd docker; make build
```

Then run the container:

```bash
docker run -v /path/to/my_data:/data -it ont-assembly-polish
```

You will be dropped into the directory `/home/ont-assembly-polish`, then simply edit `config.mk` and run the pipeline.
The `-v` flag will make the `/path/to/my_data` directory on the host available under `/data` in the container.

Application dependencies
------------------------

- [Canu](http://canu.readthedocs.io)
- [samtools](http://www.htslib.org/)
- [BWA](http://bio-bwa.sourceforge.net)
- [racon](https://github.com/isovic/racon) - the pipeline will download and build it
- [minimap](https://github.com/lh3/minimap) and [miniasm](https://github.com/lh3/miniasm)
- [pilon](https://github.com/broadinstitute/pilon/wiki) - the pipeline will download it

Evaluation on simulated data
----------------------------

In order to evaluate the performance of the pipeline we have simulated long and short reads from the yeast genome and measured the accuracy of recovered contigs
after various stages of correction.

Long reads were simulated using an in-house script under the following conditions:
- Number of reads: 150000
- Read lengths sampled from a gamma distribution with mean 6000 and shape 0.5 and a minimum read length of 600
- Simulated error rate was 0.1, errors were uncorrelated events of size one with substitution:insertion:deleltion ratio of 1:1:2

Short reads were simulated using [simLibrary and simNGS](https://www.ebi.ac.uk/goldman-srv/simNGS/):
- Simulated Illumina data consisted of paired-end reads of size 101, with the default insert length of 400
- Simulation runfile can be found under: data/s_1_4x.runfile
- The number of simulated read pairs was 21666129 (360x fragment coverage)

We have measured the accuracy of recovered contigs after various correction stages using dnadiff from the [mummer](http://mummer.sourceforge.net/) package and [last](http://last.cbrc.jp/):

![alt text](https://github.com/nanoporetech/ont-assembly-polish/blob/master/results/ddif_plots.png "dnadiff accuracies")

![alt text](https://github.com/nanoporetech/ont-assembly-polish/blob/master/results/la_plots.png "lastal accuracies")

### Conclusions from the evaluation:
- At the simulated error rate, canu alone recovers high accuracy contigs.
- Both dnadiff and lastal accuracies suggest that racon and pilon polishing increases contig accuracy.
- Lastal accuracies suggest that the effect of racon and pilon polishing is additive, though the increase in accuracy is not substantial.

License
-------

(c) 2016 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

