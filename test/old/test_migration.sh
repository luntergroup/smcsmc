echo "Migration test"
echo

echo "Building test data..."
python pop_struct.py sim-1-migration 2 1 median
echo

cmd="../smcsmc -nsam 2 -t 1000 -r 200 1000000 -I 2 1 1 -ej 0.6 1 2 -eM 0 0.1 -EM 0 -Np 1000 -seg sim-1-migrationSamples2msdata1.seg -seed 1 1 1 -p 1*3+15*4+1 -tmax 4"
echo "Running migration test..."
echo "Command line:"
echo ${cmd}
${cmd} | sed 's/..Particle[^\n]*completed.//g' > test_migration.output
rm sim-1-migrationSamples2msdata1*
diff test_migration.output test_migration.truth
