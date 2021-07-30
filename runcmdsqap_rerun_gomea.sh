export JULIA_NUM_THREADS=1
echo "Re-running QAP for GOMEAs"

cat ./cmdsqap_rerun_gomea.txt | xargs -P 6 -n 1 -d '\n' sh -c
# TODO: Add more copy commands for other problems - if needed
mkdir -p ./results/QAP && cp -r ./src/problems/QAP/results/* ./results/QAP/
tar -czf ./results.tar.gz ./results