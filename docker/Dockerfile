FROM ubuntu:16.04

LABEL version="1.0.2" \
     description="Oxford MMM CompassCompact Pipeline" \
     maintainer="crookscs.it@ndm.ox.ac.uk" \
     dockerhub="oxfordmmm/compasscompact:v1.0.2"

#Fresh update packages
RUN apt-get update && \
    apt-get install -y software-properties-common

#Update essential libraries
RUN apt-get update -y &&\
    apt-get install -y  \
    build-essential     \
    curl                \
    git                 \
    libncurses5-dev     \
    libbz2-dev          \
    liblzma-dev         \   
    nfs-common          \
    portmap             \
    python-dev          \
    python-lxml         \
    python-pip          \
    python-pycurl       \
    samtools htop       \
    sed                 \
    zip                 \
    zlib1g-dev        &&\
    apt-get clean     &&\
    apt-get purge     &&\
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#Special package for g++-4.4
RUN add-apt-repository ppa:ubuntu-toolchain-r/test &&\
    apt-get update &&\
    apt-get install -y g++-4.4 &&\
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# #Install Oracle Java 1.7
WORKDIR /opt
COPY tools/jre-7u80-linux-x64.tar.gz .
RUN  tar -xvzf jre-7u80-linux-x64.tar.gz -C . &&\
     rm "jre-7u80-linux-x64.tar.gz" &&\
     mv "jre1.7.0_80" java17 &&\
     update-alternatives --install /usr/bin/java java /opt/java17/bin/java 100
    
#Install python and libraries
RUN pip install --upgrade pip 
RUN pip install numpy biopython pysam==0.8.3

#Create compass user
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse compass && \
    echo `echo "compass\ncompass\n" | passwd password` && \
    chown compass:compass /home/compass

#Set up directories
RUN mkdir /home/compass/PIPELINE /data &&\
    chown compass:compass /data

#Define env variables and copy source code
ENV HOME=/home/compass
ENV PIPELINE=$HOME/PIPELINE
ENV PATH=${PATH}:$PIPELINE/compass
ENV COMPASSCFG=${PIPELINE}/lib/compass.cfg
WORKDIR $PIPELINE

#Compass legacy code
COPY compass compass

#Compass legacy tools
COPY gorm-tools gorm-tools

#Compile stampy
RUN cd gorm-tools/stampy-1.0.23 &&\
    sed -i 's/-Wl//' makefile &&\
    make

#Switch User
USER compass
ENTRYPOINT [ "/bin/bash" ]