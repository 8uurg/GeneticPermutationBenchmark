from itertools import product

start_indices = [1]
instance_ids = range(1, 80)
approach_ids = [1, 2]
num_exp = 20

def generate_run_cmd(start_idx, instance_id, approach_id, num_exp):
    return f"cd ./src/problems/PFS && julia -O3 --project=../../.. ./run_single_taillard_per_cat_eval_limit_focusrrg.jl {start_idx} {instance_id} {approach_id} {num_exp}\n"

with open("./cmdspfs_rerun_all.txt", "w") as f:
    for (start_idx, instance_id, approach_id) in product(start_indices, instance_ids, approach_ids):
        f.write(generate_run_cmd(start_idx, instance_id, approach_id, num_exp))