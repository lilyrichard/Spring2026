#!/bin/bash

SCRIPT1="/home/lricha95/LIED/mpcalc_TI_rate.py"
SCRIPT2="/home/lricha95/LIED/mpcalc_MF_DCS.py"

for E in $(seq 50 2 80)
do
    echo "Calculating TI Rate energy E= $E eV"
    python3 "$SCRIPT1" --E "$E" 
    echo "Running scattering energy E = $E eV"
    python3 "$SCRIPT2" --E "$E"
done

