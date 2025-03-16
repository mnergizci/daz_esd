#!/bin/bash

## Set daz_esd Path robustly
export daz_esd="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)"

## Set paths explicitly and clearly
export PATH="$daz_esd/bin:$PATH"
export PYTHONPATH="$daz_esd/lib:$PYTHONPATH"