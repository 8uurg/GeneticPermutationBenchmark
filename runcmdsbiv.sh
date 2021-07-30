export JULIA_NUM_THREADS=1

echo "Running BIV incrementally for all"

# Part of the re-run on due to
# GOMEA due to inconsistency (implementation is steady state, while it shouldn't be)
# However, as this will be ran on a different machine,
# all approaches will be re-ran instead.

cat ./cmdsbiv.txt | xargs -P 1 -n 1 -d '\n' sh -c
# TODO: Add more copy commands for other problems - if needed
mkdir -p ./results/BIV && cp -r ./src/problems/BIV/results/* ./results/BIV/
tar -czf ./results.tar.gz ./results