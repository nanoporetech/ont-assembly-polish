FROM ubuntu:18.04
MAINTAINER Botond.Sipos@nanoporetech.com

# Upgrade and install necessary packages:
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get upgrade -y &&\
    DEBIAN_FRONTEND=noninteractive apt-get install -y sudo git make build-essential mummer last-align python-numpy python-matplotlib time bwa samtools software-properties-common gnuplot\
    zlib1g-dev mc wget libatlas-base-dev python-pip python-pandas cmake

# Install latest Java:
RUN add-apt-repository ppa:webupd8team/java &&\
    apt-get update &&\
    cd /var/lib/dpkg/info &&\
    sudo sed -i 's|JAVA_VERSION=8u201|JAVA_VERSION=8u211|' oracle-java8-installer.* &&\
    sudo sed -i 's|PARTNER_URL=http://download.oracle.com/otn-pub/java/jdk/8u201-b09/42970487e3af4f5aa5bca3f542482c60/|PARTNER_URL=https://download.oracle.com/otn/java/jdk/8u211-b12/478a62b7d4e34b78b671c754eaaf38ab/|' oracle-java8-installer.* &&\
    sudo sed -i 's|SHA256SUM_TGZ="cb700cc0ac3ddc728a567c350881ce7e25118eaf7ca97ca9705d4580c506e370"|SHA256SUM_TGZ="c0b7e45330c3f79750c89de6ee0d949ed4af946849592154874d22abc9c4668d"|' oracle-java8-installer.* &&\
    sudo sed -i 's|J_DIR=jdk1.8.0_201|J_DIR=jdk1.8.0_211|' oracle-java8-installer.* &&\
    cd &&\
    (echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections) &&\
     DEBIAN_FRONTEND=noninteractive apt-get install -y oracle-java8-installer

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
