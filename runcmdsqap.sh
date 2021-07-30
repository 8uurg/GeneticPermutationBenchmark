cat ./cmdsqap.txt | xargs -P 0 -n 1 -d '\n' sh -c
# TODO: Add more copy commands for other problems - if needed
mkdir -p ./results/QAP && cp -r ./src/problems/QAP/results/* ./results/QAP/
tar -czf ./results.tar.gz ./results