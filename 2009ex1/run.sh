#!/bin/bash

ulimit -s unlimited
for file in ./in/*
do
    ./2009ex1_exec "$file" > $(echo $file | sed 's/in/out/')
done