#!/bin/bash
# This script runs GlauconiteFilterer from within the docker container
# It serves as the ENTRYPOINT of the docker.
/root/miniconda3/bin/python GlauconiteFilterer/RunGlauconiteFilterer.py -j /UserFiles/docker_json_vars.json >> /Outputfolder/output.txt 2>> /Outputfolder/error.txt