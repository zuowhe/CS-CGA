function [dag, g_best_score, conv, iterations, analysis_stats] = HSM_analysis_process(data, N, M, MP, scoring_fn, bnet, tour, p_value, variant, ci_update_rule)
if nargin < 10 || isempty(ci_update_rule)
    ci_update_rule = 'paper';
end
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
[CI_new, count_CIchange] = update_analysis_ci_threshold(ci_update_rule, p_avg, N, M, BN_NodesNum, 0, 0, 0);
CI_init = CI_new;
SuperStructure = xor(true(BN_NodesNum), diag(true(1, BN_NodesNum)));
for row_id = 1:BN_NodesNum-1
    for col_id = row_id+1:BN_NodesNum
        if p_value(row_id, col_id) > CI_init
            SuperStructure(row_id, col_id) = false;
            SuperStructure(col_id, row_id) = false;
        end
    end
end
[pop, l_map_MI] = CS_CGA_initialize_population(SuperStructure, N2, norm_MI);
pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);
score_state = HSM_analysis_init_state(variant, BN_NodesNum);
[score, score_state] = HSM_analysis_score_population(data, ns, pop, scoring_fn, score_state);
[g_best, g_best_score] = get_best(score, pop);
Dif_BIC = 0;
l_map_MI2 = l_map_MI;
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
    [pop, ~] = EEAM_first_stage_mutation(N2, l_map_MI, pop_1, Dif_BIC, M, iter);
    pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);
    [l_map_MI2] = PCR_build_relaxed_edge_map(BN_NodesNum, p_value, norm_MI, CI_new);
    diff_map = setdiff(l_map_MI2, l_map_MI, 'rows');
    [pop, now_repeats_indices] = EEAM_handle_duplicate_individuals( ...
        pop, score_state.duplicateHistory, g_best, iter, M, false);
    if ~isempty(diff_map)
        pop = EEAM_second_stage_mutation(diff_map, pop, now_repeats_indices);
    end
    pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);
    [score, score_state] = HSM_analysis_score_population(data, ns, pop, scoring_fn, score_state);
    [~, g_best_score2] = get_best(score, pop);
    Dif_BIC = score_gap_normalize(g_best_score2, g_best_score, score);
    [CI_new, count_CIchange] = update_analysis_ci_threshold( ...
        ci_update_rule, p_avg, N, M, BN_NodesNum, Dif_BIC, CI_new, count_CIchange);
    [g_best, g_best_score, pop, score] = update_elite(g_best, g_best_score, pop, score);
    conv = update_conv(conv, g_best, g_best_score, bnet.dag, iter);
    iterations = iter;
end
analysis_stats = HSM_analysis_cache_stats(score_state);
analysis_stats.finalAlpha = CI_new;
analysis_stats.initialAlpha = CI_init;
analysis_stats.alphaUpdates = count_CIchange;
analysis_stats.searchSpaceEdges = size(l_map_MI2, 1);
dag = {g_best};
end
function [CI_new, count_CIchange] = update_analysis_ci_threshold(rule, p_avg, N, M, nodes_num, Dif_BIC, CI_new, count_CIchange)
rule = lower(strrep(strrep(strtrim(rule), '_', '-'), ' ', '-'));
switch rule
    case {'paper', 'paper-step', 'paper-formula', 'updateci-paper'}
        [CI_new, count_CIchange] = PCR_update_alpha(p_avg, N, M, nodes_num, Dif_BIC, CI_new, count_CIchange);
    otherwise
        error('HSM_analysis_process:UnknownCIUpdateRule', 'Unknown CI update rule: %s', rule);
end
end
