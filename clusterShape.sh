#!/bin/bash

for rs in {1..2}; do
    for hit in {1..2}; do
        for det in {0..4}; do
            ./clusterShape ${rs} ${hit} ${det}
        done;
    done;
done;
