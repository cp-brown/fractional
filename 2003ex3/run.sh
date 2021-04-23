#!/bin/bash

for file in ./in/*
do
    ./2003ex3_exec "$file" > $(echo $file | sed 's/in/out/')
done