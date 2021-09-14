#!/bin/bash
folder=$1
tgz=${folder::-1}.tgz
set -x
tar -zcf $tgz $folder --remove-files

