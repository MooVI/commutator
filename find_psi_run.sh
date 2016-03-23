#!/bin/bash


SHORT=p:o:u
LONG=profile,output,unsafe:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
profile=""
outfile=""
if [[ $? != 0 ]]; then
    exit 2
fi
eval set -- "$PARSED"

body="find_psi_body.py"

while true; do
    case "$1" in
        -p|--profile)
            profile="-m cProfile -o $2"
            shift 2
	    ;;
	-u|--unsafe)
            body="find_psi_body_unsafe.py"
            shift 1
	    ;;
        -o|--output)
            outFile="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

if [[ $# != 1 ]]; then
    echo "$0: A single input file is required."
    exit 4
fi

cat $1 ${body} > find_psi_temp.py
exec &> "${outFile}"
python3 -u ${profile} find_psi_temp.py 
