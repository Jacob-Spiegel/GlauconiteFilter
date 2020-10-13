#!/bin/bash

# This script is an example for running GlauconiteFilterer within a docker.
# To modify the protein, pocket dimensions, and GA paremeters... please 
# create a JSON file with the desired user variables.
# An example JSON is provided at: /GlauconiteFilterer/docker/examples/sample_submit_GlauconiteFilterer_docker.json


# Make sure we are in the /GlauconiteFilterer/docker/ directory


# sudo should only be run in Linux or MacOS
# If Windows please instead just open the terminal with admin privileges
#   and omit the 'sudo'

sudo python ./glauconite_in_docker.py -j ./examples/sample_submit_GlauconiteFilterer_docker.json
