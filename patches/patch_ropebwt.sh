#!/bin/sh

fp=$1

sed -i 's/fprintf(stderr, "\[M::%s\] Turn off/; \/\/ fprintf(stderr, "\[M::%s\] Turn off/g' $fp