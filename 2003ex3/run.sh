#!/bin/bash

for file in ./in/*
do
    ./2003ex3 "$file" > $(echo $file | sed 's/in/out/')
done