#!/bin/sh


time ./dxter.x 5 | tee tensorTest/tensorResults_MP3
time ./dxter.x 6 | tee tensorTest/tensorResults_W
time ./dxter.x 7 | tee tensorTest/tensorResults_X
time ./dxter.x 8 | tee tensorTest/tensorResults_U
time ./dxter.x 9 | tee tensorTest/tensorResults_Q
time ./dxter.x 10 | tee tensorTest/tensorResults_P
time ./dxter.x 11 | tee tensorTest/tensorResults_H
time ./dxter.x 12 | tee tensorTest/tensorResults_F
time ./dxter.x 13 | tee tensorTest/tensorResults_G
time ./dxter.x 14 | tee tensorTest/tensorResults_z_small
time ./dxter.x 15 | tee tensorTest/tensorResults_Z_big
time ./dxter.x 16 | tee tensorTest/tensorResults_Tau