#!/bin/bash
WORKDIR=$PWD
cd /share/hormozdiarilab/Codes/Stella
source venv3/bin/activate
cd experiments/simulator
python -m simulator "$@"
