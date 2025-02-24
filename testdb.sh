#!/usr/bin/env bash

for f in "$@"
do 
  echo "Compreesing $f size matrix 348945080"
  ./out 348945080 $f
done
