source ../../rescomp.config.sh
python bottleneck_bylag_falsestart_experiment.py -j 10 --force -c -C "-P lunter.prjb -q long.qb -pe shmem 4"
