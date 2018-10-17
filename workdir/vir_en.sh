#!/bin/bash

#= Give execute permission to your script
#= chmod +x yourscript.sh
#= ./yourscript.sh


pip install --user virtualenv

#= create a virtual environment (venv) with name as py2env
virtualenv -p /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/python py2env

#==== activate venv
source py2env/bin/activate
#==== inside virtual environment
pip install requests
pip install -r requirements.txt

#deactivate
#==== outside virtual environment
