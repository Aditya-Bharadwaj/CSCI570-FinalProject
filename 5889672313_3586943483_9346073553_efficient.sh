#!/bin/bash

if [ -n "$1" ]; then
  echo "Running the program on $1"
  python3 5889672313_3586943483_9346073553_efficient.py $1
else
  echo "No input file provided, running on input.txt"
  python3 5889672313_3586943483_9346073553_efficient.py input.txt
fi
