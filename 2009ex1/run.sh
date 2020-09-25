#!/bin/bash

for file in ./in/*
do
    mpiexec -np 4 2009ex1 "$file" > $(echo $file | sed 's/in/out/')
done