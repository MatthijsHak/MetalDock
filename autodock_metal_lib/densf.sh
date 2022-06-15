#!/bin/bash

$AMSBIN/densf << eor

ADFFILE adf.rkf

Grid Inline


End

Potential coul scf
eor
