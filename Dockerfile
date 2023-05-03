FROM ubuntu:16.04

MAINTAINER MIKANG SIM

RUN apt-get update && apt-get install -y \
	perl \
	python3 \
	git \
	gcc \
	g++ \
	cpanminus \
	build-essential \
	pkg-config \
	libgd-dev \
	libncurses-dev \
	libghc-bzlib-dev \
	libboost-all-dev \
	build-essential \
	libz-dev \
	openjdk-8-jdk \
	openjdk-8-jre \
	make \
	zip \
	gzip \
	wget \
	vim

RUN apt-get update -y && \
	apt-get upgrade -y && \
	apt-get dist-upgrade -y && \
	apt-get install build-essential software-properties-common -y && \
	add-apt-repository ppa:ubuntu-toolchain-r/test -y && \
	apt-get update -y && \
	apt-get install gcc-7 g++-7 -y && \
	update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7 && \
	update-alternatives --config gcc


ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH

RUN wget --quiet --no-check-certificate https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh && \
    echo "45c851b7497cc14d5ca060064394569f724b67d9b5f98a926ed49b834a6bb73a *Anaconda3-2019.03-Linux-x86_64.sh" | sha256sum -c - && \
    /bin/bash /Anaconda3-2019.03-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Anaconda3-2019.03-Linux-x86_64.sh && \
    echo export PATH=$CONDA_DIR/bin:'$PATH' > /etc/profile.d/conda.sh

# Whole-genome alignment
RUN conda config --append channels conda-forge
RUN conda install -c bioconda ucsc-twobittofa
RUN conda install -y -c bioconda repeatmasker
RUN conda install -c bioconda ucsc-twobitmask
RUN conda install -c bioconda ucsc-fatotwobit
RUN conda install -c bioconda lastz
RUN conda install -c bioconda ucsc-fasize
RUN conda install -c bioconda ucsc-fasplit
RUN conda install -c bioconda ucsc-chainmergesort
RUN conda install -c bioconda ucsc-chainnet
RUN conda install -c bioconda ucsc-lavtopsl
RUN conda install -c bioconda ucsc-axtchain
RUN conda install -c bioconda ucsc-chainantirepeat
RUN conda install -c bioconda ucsc-netsyntenic
RUN conda install -c bioconda samtools
RUN conda install -c bioconda ucsc-fafilter
# Assembly
RUN conda install -c bioconda soapdenovo2
RUN conda install -c bioconda masurca
RUN conda install -c bioconda flye
# Read mapping
RUN conda install -c bioconda bwa
RUN conda install -c bioconda bowtie2
RUN conda install -c bioconda minimap2

RUN conda update --all

RUN cpanm File::Basename
RUN cpanm Parallel::ForkManager
RUN cpanm Getopt::Long
RUN cpanm FindBin
RUN cpanm List::Util
RUN cpanm POSIX
RUN cpanm Switch
RUN cpanm Role::Tiny
RUN cpanm Sub::Quote
RUN cpanm XML::SAX
RUN cpanm IO::String

RUN cpanm XML::DOM
RUN cpanm XML::LibXML
RUN cpanm XML::Parser::PerlSAX
RUN cpanm XML::Twig
RUN cpanm Test::RequiresInternet
RUN cpanm Graph::Directed
RUN cpanm Test::Most
RUN cpanm Test::Memory::Cycle
RUN cpanm YAML
RUN cpanm LWP::UserAgent
RUN cpanm XML::LibXML::Reader

RUN cpanm Bio::TreeIO
#RUN cpanm Data::Grove
#RUN cpanm LWP
#RUN cpanm Test::Simple
#RUN cpanm ExtUtils::ParseXS
#RUN cpanm Scope::Guard
#RUN cpanm Alien
#RUN cpanm HTTP::Negotiate
#RUN cpanm Module::Pluggable
#RUN cpanm Alien::Build
#RUN cpanm XML::SAX
RUN perl -MCPAN -e 'install Parallel::ForkManager'


RUN apt-get update && apt-get install -y libcurl4-openssl-dev

COPY assess /ASSESS
ENV PATH /ASSESS/src:$PATH
ENV PATH /ASSESS/ASSESS.b12:$PATH
ENV PATH /ASSESS/ASSESS.b12/code/makeBlocks:$PATH
WORKDIR /ASSESS/ASSESS.b12/code/makeBlocks
RUN make
WORKDIR /ASSESS/bin
RUN git clone https://github.com/marbl/canu.git
WORKDIR canu/src
RUN make -j 1
WORKDIR /ASSESS/bin
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
RUN tar xvfj samtools-1.16.1.tar.bz2
RUN rm samtools-1.16.1.tar.bz2
WORKDIR samtools-1.16.1
RUN ./configure --prefix=/ASSESS/bin/samtools-1.16.1
RUN make
RUN make install
WORKDIR /ASSESS/bin
RUN wget -O longranger-2.2.2.tar.gz "https://cf.10xgenomics.com/releases/genome/longranger-2.2.2.tar.gz?Expires=1664383316&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvZ2Vub21lL2xvbmdyYW5nZXItMi4yLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjY0MzgzMzE2fX19XX0_&Signature=aXe0VpDJv1yB9f6glUKQzFl~6raJS1yQuJ4uO23ho6NWTxtayCV8jdpAUbl9nVgs665ORkrlw3fu9C4miuVp288S9~zozAEAi8UGfUtO6ZkAIG0ehEOz0rndL1X9FlX8YQZuIlbNpOQpMVwvK18mrVfsIY0t43rTIKInKwamfaPPq6ME55bjqW1QVAUYrRj7jMTr08sDYzxGjGj5zlLn2kiGMFrgmJqQftts23fnqisbjDRnP-FaSZYGHYZKJ74qQDMpzxgx94d~VuZ8qACyts~6YvpfMaMBZC2CpYWBzXj~3nKcFWg2xITUcCZH9jA0nT0BWQUeG~A8zumP79kG-A__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
RUN tar -xzvf longranger-2.2.2.tar.gz
ENV PATH /ASSESS:$PATH
ENV PATH /ASSESS/bin:$PATH
ENV PATH /ASSESS/bin/longranger-2.2.2:$PATH
RUN cp /ASSESS/sources/.bashrc ~/
RUN mkdir /assess_wd
WORKDIR /assess_wd
