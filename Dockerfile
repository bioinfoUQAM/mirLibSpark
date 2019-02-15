# master progress of dockerfile
# date: 2019-02-11
# author: Chao-Jung Wu
# memo run  : docker run -it docker_pyspark3_new2 /bin/bash
# memo build: docker build --tag=docker_pyspark3_new2 .

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
#RUN \
#    pip install --trusted-host pypi.python.org -r requirements.txt && \
#    pip install pyspark --no-cache-dir


RUN \
    pip install statsmodels && \
    pip install seaborn && \
    pip install pyspark --no-cache-dir


# Install curl
RUN \
    apt-get update && \
    apt-get -y install curl


# Install RNAfold (I think C compiler comes with pyspark.)
RUN curl https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.11.tar.gz -o file.tar.gz
RUN tar xvzf file.tar.gz
RUN \ 
    cd ViennaRNA-2.4.11 && \
    ./configure && \
    make && \
    make install&& \
    cd .. && \
    rm -fr ViennaRNA-2.4.11/



# Install bowtie
#= Do not use bowtie-1.2.2, it does not compile correctly.
#= do not curl the link, it dose not download correctly
# curl https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.0/bowtie-1.2-source.zip -o xxx.zip
# tar -xvf xxx.zip
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