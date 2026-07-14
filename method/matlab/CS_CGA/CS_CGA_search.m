function [dag, g_best_score, conv, iterations] = CS_CGA_search(data, N, M, MP, scoring_fn, bnet, tour, saved_filename, p_value, ci_update_rule)
if nargin < 10 || isempty(ci_update_rule)
    ci_update_rule = 'paper';
end

node_count = size(bnet.dag, 1);
ns = bnet.node_sizes;
conv = struct('f1', zeros(1, M), 'se', zeros(1, M), 'sp', zeros(1, M), 'sc', zeros(1, M));
iterations = 0;
max_realtime = get_max_realtime(node_count);
start_time = tic;
population_size = bitsll(N, 1);

[MI, norm_MI] = get_MI(data, ns);
norm_MI(isnan(norm_MI)) = 0;

p_avg = mean(p_value(:));
fprintf('Mean p-value: %2.3f\n', p_avg);
[alpha_current, alpha_updates] = update_alpha(ci_update_rule, p_avg, N, M, node_count, 0, 0, 0);
alpha_initial = alpha_current;
fprintf('Initial alpha: %9.5f     ', alpha_initial);

superstructure = build_superstructure(p_value, alpha_initial);
[pop, reliable_edge_map] = CS_CGA_initialize_population(superstructure, population_size, norm_MI);
fprintf('Initial superstructure undirected edges: %d\n', size(reliable_edge_map, 1));
pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);

saved_file = fopen(saved_filename, 'w');
cleanup_file = onCleanup(@() fclose(saved_file));

max_cache_size = 200 * node_count;
cache = cell(1, node_count);
for k = 1:node_count
    cache{k} = struct('masks', [], 'scores', []);
end

pop_history = containers.Map('KeyType', 'char', 'ValueType', 'any');
score_history = [];
number_edges = [];
history_length = 0;
flag_cache = true;
global_cache_hits = 0;
local_cache_overwrites = 0;

[score, pop_history, score_history, number_edges, hits, history_length, flag_cache, cache, overwrites] = ...
    HSM_score_population(data, population_size, ns, pop, pop_history, score_history, number_edges, history_length, scoring_fn, flag_cache, cache, max_cache_size);
global_cache_hits = global_cache_hits + hits;
local_cache_overwrites = local_cache_overwrites + overwrites;

[g_best, g_best_score] = get_best(score, pop);
score_gap = 0;
relaxed_edge_map = reliable_edge_map;

for iter = 1:M
    if toc(start_time) > max_realtime
        iterations = iter;
        [g_best_score, conv] = CS_CGA_finalize_output(g_best, data, bnet, scoring_fn, M, conv, iterations);
        break;
    end

    [norm_score, ~] = score_normalize(score, g_best_score, false);
    if ~isempty(find(norm_score, 1))
        [offspring, ~] = selection_tournament(N, population_size, pop, score, tour);
        offspring = CS_CGA_parent_set_crossover(N, offspring);
    else
        offspring = pop;
    end

    [pop, ~] = EEAM_first_stage_mutation(population_size, reliable_edge_map, offspring, score_gap, M, iter);
    pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);

    % PCR relaxes the CI boundary only when the search state calls for exploration.
    relaxed_edge_map = PCR_build_relaxed_edge_map(node_count, p_value, norm_MI, alpha_current);
    recovered_edge_map = setdiff(relaxed_edge_map, reliable_edge_map, 'rows');
    [pop, repeated_indices] = EEAM_handle_duplicate_individuals(pop, pop_history, g_best, iter, M, true);

    if ~isempty(recovered_edge_map)
        pop = EEAM_second_stage_mutation(recovered_edge_map, pop, repeated_indices);
    end
    pop = CS_CGA_repair_DAG_by_MI(pop, MI, MP);

    [score, pop_history, score_history, number_edges, hits, history_length, flag_cache, cache, overwrites] = ...
        HSM_score_population(data, population_size, ns, pop, pop_history, score_history, number_edges, history_length, scoring_fn, flag_cache, cache, max_cache_size);
    global_cache_hits = global_cache_hits + hits;
    local_cache_overwrites = local_cache_overwrites + overwrites;

    [~, generation_best_score] = get_best(score, pop);
    score_gap = score_gap_normalize(generation_best_score, g_best_score, score);
    [alpha_current, alpha_updates] = update_alpha(ci_update_rule, p_avg, N, M, node_count, score_gap, alpha_current, alpha_updates);

    [g_best, g_best_score, pop, score] = update_elite(g_best, g_best_score, pop, score);
    conv = update_conv(conv, g_best, g_best_score, bnet.dag, iter);
    fprintf(saved_file, '%d\n', g_best_score);
    iterations = iter;
end

fprintf('Final alpha: %f     ', alpha_current);
fprintf('Final superstructure undirected edges: %d\n', size(relaxed_edge_map, 1));
fprintf('HSM global cache hits: %d      ', global_cache_hits);
fprintf('HSM local cache overwrites: %d      ', local_cache_overwrites);
fprintf('PCR alpha updates: %d\n', alpha_updates);

dag = {g_best};
end

function superstructure = build_superstructure(p_value, alpha)
node_count = size(p_value, 1);
superstructure = xor(true(node_count), diag(true(1, node_count)));
for i = 1:node_count-1
    for j = i+1:node_count
        if p_value(i, j) > alpha
            superstructure(i, j) = false;
            superstructure(j, i) = false;
        end
    end
end
end

function [alpha_current, alpha_updates] = update_alpha(rule, p_avg, N, M, node_count, score_gap, alpha_current, alpha_updates)
rule = lower(strrep(strrep(strtrim(rule), '_', '-'), ' ', '-'));
switch rule
    case {'paper', 'paper-step', 'paper-formula', 'pcr'}
        [alpha_current, alpha_updates] = PCR_update_alpha(p_avg, N, M, node_count, score_gap, alpha_current, alpha_updates);
    otherwise
        error('Unknown PCR update rule: %s', rule);
end
end
