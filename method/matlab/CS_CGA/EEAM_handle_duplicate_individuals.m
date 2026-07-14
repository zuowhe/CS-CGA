function [pop, now_repeats_indices] = EEAM_handle_duplicate_individuals(pop, pop_history, g_best, iter, max_iter, enable_logging)
if nargin < 6, enable_logging = true; end
repeats_info = struct();
[~, now_repeats_indices] = CS_CGA_find_duplicate_individuals(pop);
now_repeats_indices = now_repeats_indices(:)';
repeats_info.total_current = length(now_repeats_indices);
common_repeated_indices = [];
if isempty(now_repeats_indices) || isempty(pop_history) || ~isa(pop_history, 'containers.Map')
    common_repeated_indices = [];
else
    for idx = now_repeats_indices
        key = CS_CGA_dag_key(pop{idx});
        if isKey(pop_history, key)
            common_repeated_indices(end + 1) = idx;
        end
    end
end
repeats_info.common_count = length(common_repeated_indices);
repeats_info.replaced_count = 0;
repeats_info.indices = [];
if ~isempty(common_repeated_indices)
    replace_ratio = 0.5 * (iter - 1) / (max_iter - 1);
    num_to_replace = max(1, round(length(common_repeated_indices) * replace_ratio));
    selected_idx = common_repeated_indices(randperm(length(common_repeated_indices), num_to_replace));
    for k = 1:length(selected_idx)
        pop(selected_idx(k)) = {g_best};
    end
    repeats_info.replaced_count = length(selected_idx);
    repeats_info.indices = selected_idx;
end
if enable_logging
    log_iterations = [10, 50, 100, 190];
    if ismember(iter, log_iterations)
        fprintf('Iteration %d duplicate individuals: %d; ', iter, repeats_info.total_current);
    end
end
end
