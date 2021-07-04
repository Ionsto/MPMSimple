#!/bin/bash
g++ *.cpp -O3 -ffast-math -std=c++2a -fopenmp && ./a.out && ffmpeg -framerate 60  -i './Frames/frame %d .png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4 -y
