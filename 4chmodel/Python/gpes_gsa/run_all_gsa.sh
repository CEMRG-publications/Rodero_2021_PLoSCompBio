#!/bin/bash

for i in `seq 0 7`
do
	python3 run_gsa.py marina ${i}
done