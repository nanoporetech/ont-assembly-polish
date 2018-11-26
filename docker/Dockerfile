FROM ubuntu:16.04
MAINTAINER Botond.Sipos@nanoporetech.com

# Upgrade and install necessary packages:
RUN apt-get update && apt-get upgrade -y &&\
    apt-get install -y sudo git make build-essential mummer last-align python-numpy python-matplotlib time bwa samtools software-properties-common gnuplot\
    zlib1g-dev mc wget libatlas-base-dev python-pip python-pandas

# Install latest Java:
RUN add-apt-repository ppa:webupd8team/java &&\
    apt-get update &&\
    (echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections) &&\
     apt-get install -y oracle-java8-installer

# Start working in /opt
WORKDIR /opt

# Compile and install canu:
RUN git clone https://github.com/marbl/canu.git &&\
    cd canu/src && make -j && cd -
ENV PATH=/opt/canu/Linux-amd64/bin:$PATH
    
# Compile and install minimap:
RUN git clone https://github.com/lh3/minimap2 && (cd minimap2 && make) &&\
    cp minimap2/minimap2 /usr/local/bin && rm -r minimap2

# Compile and install miniasm:
RUN git clone https://github.com/lh3/miniasm && (cd miniasm && make) &&\
    cp miniasm/miniasm /usr/local/bin/ && rm -r miniasm

# Compile and install simNGS, a tool for simulating Illumina sequencing:
RUN git clone https://github.com/timmassingham/simNGS.git && (cd simNGS/src && make -f Makefile.linux) &&\
    cp simNGS/bin/* /usr/local/bin/ && rm -r simNGS

# Change working directory to home:
WORKDIR /home

# Clone the ont-assembly-polish project:
ARG CACHEBUST
RUN DUMMY=${CACHEBUST} git clone https://github.com/nanoporetech/ont-assembly-polish.git

# Change into the project directory:
WORKDIR /home/ont-assembly-polish
