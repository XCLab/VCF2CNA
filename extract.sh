#!/bin/bash

BASE_DIR=$PWD
tar -jxvf prep.tar.bz2
cd gc
tar -jxvf GRCh37.gc.tar.bz2
cd ../mapability
tar -jxvf hg19.map.tar.bz2
