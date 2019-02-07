#!/bin/bash

awk 'BEGIN{OFS="\t"}{ print $4,$1,$2,$3,"." }' $1

