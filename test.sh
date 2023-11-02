#!/bin/sh

cc test.c -O3 -fno-math-errno -msse4
./a.out


