cat ./cmdspfs.txt | xargs -P 0 -n 1 -d '\n' sh -c
# TODO: Add more copy commands for other problems - if needed
mkdir -p ./results/PFS && cp -r ./src/problems/PFS/results/* ./results/PFS/
tar -czf ./results.tar.gz ./results