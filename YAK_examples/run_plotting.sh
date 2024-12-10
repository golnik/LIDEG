#!/bin/bash

graphene_prog_path="../graphene"

#input file
input=input.ini

#extract number of steps
Nt=`cat $input | grep "Nt" | awk '{print $NF}'`

for (( i=1; i<=$Nt; i++ )); do
#for i in {55,}; do
    echo "step $i"
    tasks=pulse,kspace,rspace    #specify plot tasks
    python3 $graphene_prog_path/python/main.py -input $input -tstep $i -tasks $tasks #-output plot.pdf
done

exit 0
