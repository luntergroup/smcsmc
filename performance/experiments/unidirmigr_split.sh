# run simulation
# --force: overwrite any results already in database (but do not overwrite existing data in data-zigzag/)
#
# Add "-c" to submit jobs to a compute cluster using qsub
# Add "-j 60" to run 4 (models) * 15 (replicates) experiments in parallel
python unidirmigr_split.py --db udmigr-split.db --datapath data-udmigr-split --force

# make plot
python unidirmigr_split.py --db udmigr-split.db --datapath data-udmigr-split --plot
