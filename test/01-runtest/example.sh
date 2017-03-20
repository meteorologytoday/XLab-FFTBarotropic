#!/bin/bash

mkdir input
mkdir output
mkdir img

makefield-elliptic-vortex.out
main.out

python plot/draw_figs.py --input-dir=output --output-dir=img
