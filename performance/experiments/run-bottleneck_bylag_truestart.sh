source ../../rescomp.config.sh
python bottleneck_bylag_truestart_experiment.py -j 100 --force -c -C "-P lunter.prjb -q long.qb -pe shmem 4"
