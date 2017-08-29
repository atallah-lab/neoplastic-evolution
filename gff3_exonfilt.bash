#!/bin/bash

grep exon $1 > tmp
cut -f1,4,5 tmp > $1.filt
rm tmp

