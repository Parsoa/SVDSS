#!/bin/bash
module load picard-tools/2.23.4
picard ValidateSamFile -I $1 --MODE VERBOSE --MAX_OUTPUT 1000000000
