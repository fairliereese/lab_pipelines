set -e 

fname=$1
runs_dir=/runs/

# unzip
jar xf ${runs_dir}${fname}

# remove zip
rm ${runs_dir}${fname} 
