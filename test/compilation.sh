#!/bin/bash

FILE_NAME="mergesort"
CPP_FILE="mergesort.cpp"

mpic++ -O3 -std=c++11 -pedantic-errors -Wall $CPP_FILE -o $FILE_NAME
