#! /bin/bash

for alpha1 in $(seq 0.03 0.001 0.045); do
        echo -n "$alpha1: "
        run_get_pitch $alpha1 | grep TOTAL
done