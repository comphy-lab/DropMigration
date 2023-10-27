#!/bin/bash

scp -r get* $1/
scp -r *.py $1/

cd $1
echo "Post processing for $1"
python ActiveStuff_VP.py --num_workers 4 --tSnap 0.1 --L0 16.0
cd ../