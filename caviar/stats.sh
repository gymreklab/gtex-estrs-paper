#!/bin/bash

cat /storage/mgymrek/gtex-estrs/revision/caviar/output/*.tab | grep -v "num" | awk '{print $NF}' | datamash mean 1
