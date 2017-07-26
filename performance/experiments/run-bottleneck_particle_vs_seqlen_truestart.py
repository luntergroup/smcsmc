source ../../rescomp.config.sh
python bottleneck_particle_vs_seqlen_truestart.py -j 100 --force -c -C "-P lunter.prjb -q short.qb -pe shmem 4"
