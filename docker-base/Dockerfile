FROM debian:stretch

RUN apt-get update -qq \
    && apt-get install -yq --no-install-recommends \
        bash-completion \
        build-essential \
        ca-certificates \
        git \
        mercurial \
        libhdf5-dev \
        libxml2-dev \
        ssh \
        wget \
        g++ \
        gfortran \
        mercurial \
        python \
        python-lxml \
        python-h5py \
        python-numpy \
        python-scipy \
        cmake \
        mpich \
        libmpich-dev \
        libhdf5-dev \
        rsync \
        vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ARG PETSC_CONFIG_OPTIONS="--prefix=/petsc-3.5.4 --with-debugging=0 --with-shared-libraries --download-mumps=1 --download-fblaslapack=1 --download-parmetis=1 --download-metis=1 --download-scalapack=1 --download-blacs=1 --download-superlu_dist=1 --with-mpi=1"

RUN wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.5.4.tar.gz \
    && tar -xvf petsc-3.5.4.tar.gz \
    && mv petsc-3.5.4 petsc-3.5.4-src \
    && cd petsc-3.5.4-src \
    && PETSC_DIR=$PWD ./config/configure.py ${PETSC_CONFIG_OPTIONS} \
    && make PETSC_DIR=/petsc-3.5.4-src PETSC_ARCH=arch-linux2-c-opt all \
    && make PETSC_DIR=/petsc-3.5.4-src PETSC_ARCH=arch-linux2-c-opt install \
    && rm -rf petsc-3.5.4.tar.gz \
    && rm -rf petsc-3.5.4-src

ARG CONFIG_OPTIONS="--hdf5-lib-dir=/usr/lib/x86_64-linux-gnu/hdf5/serial \
        --hdf5-inc-dir=/usr/include/hdf5/serial \
        --petsc-dir=/petsc-3.5.4 \
        --with-debugging=0 \
        --cflags='-Wno-error=format-security'"

ARG UW_GIT_URL="https://github.com/LukeMondy/EarthByte_Underworld.git"

RUN git clone --depth=1 ${UW_GIT_URL}  \
    && cd /EarthByte_Underworld \
    && ./configure.py ${CONFIG_OPTIONS} \
    && ./scons.py \
    && export UW_DIR=/EarthByte_Underworld/build \
    && cd /EarthByte_Underworld/earthbyte_additions \
    && ./configure.py ${CONFIG_OPTIONS}\
    && ./scons.py \
    && cd /EarthByte_Underworld/lecode_tools \
    && ./configure.py ${CONFIG_OPTIONS}\
    && ./scons.py 

RUN mkdir /models \
    && chmod o+rwx /models \
    && cd /models \
    && git clone https://github.com/LukeMondy/lithospheric_modelling_recipe.git \
    && chmod o+rwx lithospheric_modelling_recipe \
