# master progress of dockerfile
# date: 2019-02-11
# author: Chao-Jung Wu
# version: 00.01.13
# memo run  : docker run -it environmentv13 /bin/bash
# memo build: docker build --tag=environmentv13 .
#
#
# description: python2, pysaprk, without hadoop or nodes configuration

FROM ubuntu:latest


# Copy the current directory contents into the container
COPY . /

# Install python2.7
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
    build-essential \
    ca-certificates \
    gcc \
    git \
    libpq-dev \
    make \
    python-pip \
    python2.7 \
    python2.7-dev \
    ssh \
    && apt-get autoremove \
    && apt-get clean


# Install OpenJDK 8
RUN \
  apt-get install -y openjdk-8-jdk && \
  rm -rf /var/lib/apt/lists/*


# Install PySpark and requirements
RUN \
    pip install statsmodels && \
    pip install seaborn && \
    pip install pyspark --no-cache-dir


# Install curl
RUN \
    apt-get update && \
    apt-get -y install curl

# Install RNAfold
RUN \ 
    cd ViennaRNA-2.4.11 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -fr ViennaRNA-2.4.11/


# Install bowtie
RUN \
    cd bowtie-1.2 && \
    make NO_TBB=1 && \
    make install && \
    cd .. && \
    rm -fr bowtie-1.2


COPY dustmasker /usr/local/bin/
# test: echo $'>seq1\nCGTGGCTATGATAGCGATATTCGTTTTTTT' | dustmasker


## Define working directory
#WORKDIR /data

# Define default command
CMD ["bash"]