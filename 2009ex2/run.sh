#!/bin/bash

for file in ./in/*
do
    ./2009ex2 "$file" > $(echo $file | sed 's/in/out/')
done