function [score, state, batch_stats] = HSM_analysis_score_population(data, ns, population, scoring_fn, state)
population_size = numel(population);
score = zeros(1, population_size);
batch_stats = empty_batch_stats();
for individual_id = 1:population_size
    dag = logical(population{individual_id});
    dag_key = CS_CGA_dag_key(dag);
    switch lower(state.mode)
        case 'bnt'
            [score(individual_id), state, local_stats] = score_one_with_bnt_cache(data, ns, dag, scoring_fn, state);
            batch_stats = add_local_stats(batch_stats, local_stats);
        case 'gb'
            state.globalQueries = state.globalQueries + 1;
            batch_stats.globalQueries = batch_stats.globalQueries + 1;
            if isKey(state.globalCache, dag_key)
                score(individual_id) = state.globalCache(dag_key);
                state.globalHits = state.globalHits + 1;
                batch_stats.globalHits = batch_stats.globalHits + 1;
            else
                [score(individual_id), state, local_stats] = score_one_with_bnt_cache(data, ns, dag, scoring_fn, state);
                state.globalCache(dag_key) = score(individual_id);
                state.globalKeyChars = state.globalKeyChars + length(dag_key);
                batch_stats = add_local_stats(batch_stats, local_stats);
            end
        case 'large-hash'
            state.globalQueries = state.globalQueries + 1;
            batch_stats.globalQueries = batch_stats.globalQueries + 1;
            if isKey(state.globalCache, dag_key)
                score(individual_id) = state.globalCache(dag_key);
                state.globalHits = state.globalHits + 1;
                batch_stats.globalHits = batch_stats.globalHits + 1;
            else
                [score(individual_id), state, local_stats] = score_one_with_large_hash(data, ns, dag, scoring_fn, state);
                state.globalCache(dag_key) = score(individual_id);
                state.globalKeyChars = state.globalKeyChars + length(dag_key);
                batch_stats = add_local_stats(batch_stats, local_stats);
            end
        case 'hsm'
            state.globalQueries = state.globalQueries + 1;
            batch_stats.globalQueries = batch_stats.globalQueries + 1;
            if isKey(state.globalCache, dag_key)
                score(individual_id) = state.globalCache(dag_key);
                state.globalHits = state.globalHits + 1;
                batch_stats.globalHits = batch_stats.globalHits + 1;
            else
                [score(individual_id), state, local_stats] = score_one_with_hsm_cache(data, ns, dag, scoring_fn, state);
                state.globalCache(dag_key) = score(individual_id);
                state.globalKeyChars = state.globalKeyChars + length(dag_key);
                batch_stats = add_local_stats(batch_stats, local_stats);
            end
        otherwise
            error('HSM_analysis_score_population:UnknownMode', 'Unknown cache mode: %s', state.mode);
    end
    state.duplicateHistory(dag_key) = true;
end
end
function batch_stats = empty_batch_stats()
batch_stats = struct();
batch_stats.globalQueries = 0;
batch_stats.globalHits = 0;
batch_stats.localQueries = 0;
batch_stats.localHits = 0;
batch_stats.overwrites = 0;
end
function batch_stats = add_local_stats(batch_stats, local_stats)
batch_stats.localQueries = batch_stats.localQueries + local_stats.localQueries;
batch_stats.localHits = batch_stats.localHits + local_stats.localHits;
batch_stats.overwrites = batch_stats.overwrites + local_stats.overwrites;
end
function [score, state, local_stats] = score_one_with_bnt_cache(data, ns, dag, scoring_fn, state)
node_count = size(data, 1);
node_type = repmat({'tabular'}, 1, node_count);
params = cell(1, node_count);
for node_id = 1:node_count
    params{node_id} = {'prior_type', 'dirichlet', 'dirichlet_weight', 1};
end
discrete = 1:node_count;
score = 0;
local_stats = empty_batch_stats();
for node_id = 1:node_count
    ps = parents(dag, node_id);
    local_stats.localQueries = local_stats.localQueries + 1;
    state.localQueries = state.localQueries + 1;
    [found, ~] = score_find_in_cache(state.bntCache, node_id, ps, scoring_fn);
    if found
        local_stats.localHits = local_stats.localHits + 1;
        state.localHits = state.localHits + 1;
    end
    [family_score, state.bntCache] = score_family( ...
        node_id, ps, node_type{node_id}, scoring_fn, ns, discrete, data, params{node_id}, state.bntCache);
    score = score + family_score;
end
end
function [score, state, local_stats] = score_one_with_large_hash(data, ns, dag, scoring_fn, state)
node_count = size(data, 1);
node_type = repmat({'tabular'}, 1, node_count);
params = cell(1, node_count);
for node_id = 1:node_count
    params{node_id} = {'prior_type', 'dirichlet', 'dirichlet_weight', 1};
end
discrete = 1:node_count;
score = 0;
local_stats = empty_batch_stats();
for node_id = 1:node_count
    ps = parents(dag, node_id);
    local_key = local_parent_key(node_id, ps, scoring_fn);
    local_stats.localQueries = local_stats.localQueries + 1;
    state.localQueries = state.localQueries + 1;
    if isKey(state.localHash, local_key)
        family_score = state.localHash(local_key);
        local_stats.localHits = local_stats.localHits + 1;
        state.localHits = state.localHits + 1;
    else
        [family_score, ~] = score_family( ...
            node_id, ps, node_type{node_id}, scoring_fn, ns, discrete, data, params{node_id}, []);
        state.localHash(local_key) = family_score;
        state.localKeyChars = state.localKeyChars + length(local_key);
    end
    score = score + family_score;
end
end
function [score, state, local_stats] = score_one_with_hsm_cache(data, ns, dag, scoring_fn, state)
node_count = size(data, 1);
node_type = repmat({'tabular'}, 1, node_count);
params = cell(1, node_count);
for node_id = 1:node_count
    params{node_id} = {'prior_type', 'dirichlet', 'dirichlet_weight', 1};
end
discrete = 1:node_count;
score = 0;
local_stats = empty_batch_stats();
for node_id = 1:node_count
    ps = parents(dag, node_id);
    local_stats.localQueries = local_stats.localQueries + 1;
    state.localQueries = state.localQueries + 1;
    [found, family_score] = find_hsm_parent_cache(state.localCache{node_id}, ps, node_count);
    if found
        local_stats.localHits = local_stats.localHits + 1;
        state.localHits = state.localHits + 1;
    else
        before_overwrites = get_node_overwrites(state.localCache{node_id});
        [family_score, state.localCache, ~, ~] = HSM_score_family( ...
            node_id, ps, node_type{node_id}, ns, discrete, data, params{node_id}, ...
            state.localCache, state.parentCacheSize);
        after_overwrites = get_node_overwrites(state.localCache{node_id});
        overwrite_delta = after_overwrites - before_overwrites;
        local_stats.overwrites = local_stats.overwrites + overwrite_delta;
        state.overwrites = state.overwrites + overwrite_delta;
    end
    score = score + family_score;
end
end
function [found, score] = find_hsm_parent_cache(node_cache, ps, node_count)
found = false;
score = 0;
if isempty(node_cache) || ~isfield(node_cache, 'masks') || isempty(node_cache.masks)
    return
end
mask = false(1, node_count);
if ~isempty(ps)
    mask(unique(ps)) = true;
end
for row_id = 1:size(node_cache.masks, 1)
    if isequal(node_cache.masks(row_id, :), mask)
        found = true;
        score = node_cache.scores(row_id);
        return
    end
end
end
function count = get_node_overwrites(node_cache)
count = 0;
if ~isempty(node_cache) && isfield(node_cache, 'overwrite_count')
    count = node_cache.overwrite_count;
end
end
function key = local_parent_key(node_id, ps, scoring_fn)
ps = sort(unique(ps));
if isempty(ps)
    parent_part = 'none';
else
    parent_part = sprintf('%d_', ps);
end
key = sprintf('%d|%s|%s', node_id, scoring_fn, parent_part);
end
