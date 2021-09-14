for i in `ls -d 19050*/`; do tar -zcf ${i::-1}.tgz $i --remove-files; done;
