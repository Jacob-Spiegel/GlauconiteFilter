# Get the base image from ubuntu
FROM ubuntu:18.04

MAINTAINER Jacob Spiegel <jspiegel@pitt.edu>

# Install programs through apt-get
RUN apt-get -yqq update
RUN apt-get -fyqq install
RUN apt-get -yqq install wget
RUN apt-get -yqq update
RUN apt-get -yqq install zip
# RUN apt-get -yqq install git

# Install python
RUN apt-get -yqq install bzip2
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod -R a+rwx Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b
RUN rm ./Miniconda3-latest-Linux-x86_64.sh
RUN echo "alias python=/root/miniconda3/bin/python" >> /root/.bashrc
RUN echo "alias conda=/root/miniconda3/bin/conda" >> /root/.bashrc
RUN echo "alias pip=/root/miniconda3/bin/pip" >> /root/.bashrc


# Get python dependencies
RUN /root/miniconda3/bin/conda install -y -c conda-forge rdkit rdkit=2020.03.1
RUN /root/miniconda3/bin/conda install -y numpy numpy=1.18.1
RUN /root/miniconda3/bin/conda install -y scipy scipy=1.4.1
RUN /root/miniconda3/bin/pip install func_timeout==4.3.5

# Copy over GlauconiteFilterer files
RUN mkdir /GlauconiteFilterer/
ADD ./GlauconiteFilterer /GlauconiteFilterer

# Copy over User files for GlauconiteFilterer Run
RUN mkdir /UserFiles
ADD ./temp_user_files /UserFiles

# Copy over the GlauconiteFilterer run script
ADD ./run_glauconite_in_container.bash /glauconite/run_glauconite_in_container.bash
RUN chmod -R a+rwx /glauconite/run_glauconite_in_container.bash
ADD ./run_glauconite_in_container_windows.bash /glauconite/run_glauconite_in_container_windows.bash
RUN chmod -R a+rwx /glauconite/run_glauconite_in_container_windows.bash

# Make directories for mounted host file systems
RUN mkdir /Outputfolder


RUN ls GlauconiteFilterer
ENTRYPOINT ["bash", "/glauconite/run_glauconite_in_container.bash"]

# Windows version. Automatically turned on by glauconite_in_docker.py
# ENTRYPOINT ["bash", "/glauconite/run_glauconite_in_container.bash"]


# ENTRYPOINT ["bash"]