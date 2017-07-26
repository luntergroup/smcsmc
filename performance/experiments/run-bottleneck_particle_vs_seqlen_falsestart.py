source ../../rescomp.config.sh
python bottleneck_particle_vs_seqlen_falsestart.py -j 100 --force -c -C "-P lunter.prjb -q long.qb -pe shmem 4"
