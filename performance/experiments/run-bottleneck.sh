source ../../rescomp.config.sh
python bottleneck_inference_guiding.py -j 100 --force -c -C "-P lunter.prjb -q long.qb -pe shmem 4"
