function [dag, g_best_score, conv, iterations] = CS_CGA_ablation_search(data, N, M, MP, scoring_fn, bnet, tour, saved_filename, p_value, variant)
cfg = CS_CGA_ablation_config(variant);
BN_NodesNum = size(bnet.dag, 1);
ns = bnet.node_sizes;
conv = struct('f1', zeros(1, M), 'se', zeros(1, M), 'sp', zeros(1, M), 'sc', zeros(1, M));
iterations = 0;
max_realtime = get_max_realtime(BN_NodesNum);
inStart = tic;
N2 = bitsll(N, 1);
[MI, norm_MI] = get_MI(data, ns);
norm_MI(isnan(norm_MI)) = 0;
p_avg = mean(p_value(:));
count_CIchange = 0;
if cfg.use_pcr
    [CI_new, count_CIchange] = PCR_update_alpha(p_avg, N, M, BN_NodesNum, 0, 0, 0);
else
    CI_new = cfg.fixed_alpha;
end
SuperStructure = build_superstructure(p_value, CI_new);
if isempty(find(SuperStructure, 1))
    empty_dag = false(BN_NodesNum);
    [g_best_score, conv] = CS_CGA_finalize_output(empty_dag, data, bnet, scoring_fn, M, conv, 1);
    dag = {empty_dag};
    return
end
[pop, l_map_base] = CS_CGA_initialize_population(SuperStructure, N2, norm_MI);
pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);
saved_file = fopen(saved_filename, 'w');
cleanup_file = onCleanup(@() fclose_if_open(saved_file));
max_cache_size = 200 * BN_NodesNum;
cache = cell(1, BN_NodesNum);
for k = 1:BN_NodesNum
    cache{k} = struct('masks', [], 'scores', []);
end
pop_history = containers.Map('KeyType', 'char', 'ValueType', 'any');
score_history = [];
number_edges = [];
length_history = 0;
flag_cache = true;
repeated_count = 0;
total_overwrite = 0;
[score, pop_history, score_history, number_edges, repeated_count1, length_history, flag_cache, cache, total_once_overwrite] = ...
    HSM_score_population(data, N2, ns, pop, pop_history, score_history, number_edges, length_history, scoring_fn, flag_cache, cache, max_cache_size);
repeated_count = repeated_count + repeated_count1;
total_overwrite = total_overwrite + total_once_overwrite;
[g_best, g_best_score] = get_best(score, pop);
Dif_BIC = 0;
for iter = 1:M
    if toc(inStart) > max_realtime
        iterations = iter;
        [g_best_score, conv] = CS_CGA_finalize_output(g_best, data, bnet, scoring_fn, M, conv, iterations);
        break;
    end
    [norm_score, ~] = score_normalize(score, g_best_score, false);
    if ~isempty(find(norm_score, 1))
        [pop_1, ~] = selection_tournament(N, N2, pop, score, tour);
        pop_1 = CS_CGA_parent_set_crossover(N, pop_1);
    else
        pop_1 = pop;
    end
    if cfg.use_pcr
        l_map_current = PCR_build_relaxed_edge_map(BN_NodesNum, p_value, norm_MI, CI_new);
    else
        l_map_current = l_map_base;
    end
    if cfg.use_second_stage
        stage1_map = l_map_base;
    else
        stage1_map = l_map_current;
    end
    if cfg.use_mi_mutation
        [pop, ~] = EEAM_first_stage_mutation(N2, stage1_map, pop_1, Dif_BIC, M, iter);
    else
        pop = CS_CGA_random_bitflip_mutation(N2, stage1_map, pop_1);
    end
    pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);
    if cfg.use_second_stage
        diff_map = setdiff(l_map_current, l_map_base, 'rows');
        [pop, repeat_indices] = EEAM_handle_duplicate_individuals(pop, pop_history, g_best, iter, M, false);
        if ~isempty(diff_map) && ~isempty(repeat_indices)
            if cfg.use_second_stage_mi
                pop = EEAM_second_stage_mutation(diff_map, pop, repeat_indices);
            else
                pop = random_recovered_mutation(diff_map, pop, repeat_indices);
            end
            pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);
        end
    end
    [score, pop_history, score_history, number_edges, repeated_count1, length_history, flag_cache, cache, total_once_overwrite] = ...
        HSM_score_population(data, N2, ns, pop, pop_history, score_history, number_edges, length_history, scoring_fn, flag_cache, cache, max_cache_size);
    repeated_count = repeated_count + repeated_count1;
    total_overwrite = total_overwrite + total_once_overwrite;
    [~, g_best_score2] = get_best(score, pop);
    Dif_BIC = score_gap_normalize(g_best_score2, g_best_score, score);
    if cfg.use_pcr
        [CI_new, count_CIchange] = PCR_update_alpha(p_avg, N, M, BN_NodesNum, Dif_BIC, CI_new, count_CIchange);
    end
    [g_best, g_best_score, pop, score] = update_elite(g_best, g_best_score, pop, score);
    conv = update_conv(conv, g_best, g_best_score, bnet.dag, iter);
    fprintf(saved_file, '%d\n', g_best_score);
    iterations = iter;
end
fprintf('Ablation %s final alpha: %f; repeated score hits: %d; local overwrites: %d; alpha updates: %d\n', ...
    cfg.output_name, CI_new, repeated_count, total_overwrite, count_CIchange);
dag = {g_best};
end
function ss = build_superstructure(p_value, alpha)
n = size(p_value, 1);
ss = xor(true(n), diag(true(1, n)));
for i = 1:n-1
    for j = i+1:n
        if p_value(i, j) > alpha
            ss(i, j) = false;
            ss(j, i) = false;
        end
    end
end
end
function p = random_recovered_mutation(l_map, p, indices)
l_cnt = size(l_map, 1);
if l_cnt == 0 || isempty(indices)
    return
end
m = 1 / l_cnt;
for l = 1:l_cnt
    j = l_map(l, 1);
    k = l_map(l, 2);
    for idx = indices(:)'
        l_val = get_allele(p{idx}(j, k), p{idx}(k, j));
        if m >= rand
            l_val_new = mod(l_val + round(rand), 3) + 1;
            switch l_val_new
                case 1
                    p{idx}(j, k) = false;
                    p{idx}(k, j) = false;
                case 2
                    p{idx}(j, k) = false;
                    p{idx}(k, j) = true;
                case 3
                    p{idx}(j, k) = true;
                    p{idx}(k, j) = false;
            end
        end
    end
end
end
function fclose_if_open(fid)
if isnumeric(fid) && fid > 0
    fclose(fid);
end
end
