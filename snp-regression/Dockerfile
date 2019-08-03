FROM ubuntu:16.04

# Update necessary packages
RUN apt-get update && apt-get install -qqy \
    awscli \
    python3 \
    python3-pip \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

# Install tabix
RUN mkdir /sources
WORKDIR /sources
RUN wget -O htslib-1.8.tar.bz2 https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
RUN tar -xjvf htslib-1.8.tar.bz2
WORKDIR htslib-1.8/
RUN ./configure && make && make install

RUN pip3 install numpy pandas scipy==1.2.1 statsmodels pytabix

# Add linear regression script
RUN mkdir /scripts
RUN mkdir /scratch
ADD LinRegAssociationTest_v2.py /scripts/
ADD run.sh /scripts/
WORKDIR /scripts
ENTRYPOINT ["/scripts/run.sh"]