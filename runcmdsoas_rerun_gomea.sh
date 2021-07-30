export JULIA_NUM_THREADS=10
# Rerunning GOMEA due to inconsistency (implementation is steady state, while it shouldn't be)

cat ./cmdsoas_rerun_gomea.txt | xargs -P 1 -n 1 -d '\n' sh -c
# TODO: Add more copy commands for other problems - if needed
mkdir -p ./results/OAS && cp -r ./src/problems/OAS/results/* ./results/OAS/
tar -czf ./results.tar.gz ./results