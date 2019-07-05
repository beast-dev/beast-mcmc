#!/usr/bin/env bash

function traverse(){
    cd $1
    for file in ./*Test.java
    do
        echo file #TODO: check if junit is runnable then run
    done

    for d in ./
    do
        traverse d
    done
}

traverse ../src/test.dr
