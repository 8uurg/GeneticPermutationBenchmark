export JULIA_NUM_THREADS=1
# Rerunning GOMEA due to inconsistency (implementation is steady state, while it shouldn't be)

echo "Re-running Permutation Flowshop for All"

cat ./cmdspfs.txt | xargs -P 32 -n 1 -d '\n' sh -c
# TODO: Add more copy commands for other problems - if needed
mkdir -p ./results/PFS && cp -r ./src/problems/PFS/results/* ./results/PFS/
tar -czf ./results.tar.gz ./results