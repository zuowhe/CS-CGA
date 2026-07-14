function [score, cache, cache_hit_count, overwrite_count] = HSM_score_family(j, ps, node_type, ns, discrete, data, args, cache, max_cache_size)
if nargin < 10 || isempty(max_cache_size)
    max_cache_size = 100;
end
if nargin < 9 || isempty(cache)
    n_nodes = size(data, 1);
    cache = cell(1, n_nodes);
    for k = 1:n_nodes
        cache{k} = struct('masks', [], 'scores', []);
    end
end
cache_hit_count = 0;
overwrite_count = 0;
misv = -9999;
ccc = iscell(data);
if ccc
    data = bnt_to_mat(data, misv);
end
ps = unique(ps);
[n, ncases] = size(data);
[found, score] = check_cache(cache{j}, ps);
if found
    cache_hit_count = 1;
    return;
end
dag = zeros(n);
if ~isempty(ps)
    dag(ps, j) = 1;
    ps = sort(ps);
end
bnet = mk_bnet(dag, ns, 'discrete', discrete);
fname = sprintf('%s_CPD', node_type);
if isempty(args)
    bnet.CPD{j} = feval(fname, bnet, j);
else
    bnet.CPD{j} = feval(fname, bnet, j, args{:});
end
fam = [ps j];
[tmp, available_case] = find(data(fam,:) == misv);
available_case = setdiff(1:ncases, available_case);
bnet.CPD{j} = learn_params(bnet.CPD{j}, fam, data(:, available_case), ns, bnet.cnodes);
L = log_prob_node(bnet.CPD{j}, data(j, available_case), data(ps, available_case));
S = struct(bnet.CPD{j});
score = L - 0.5 * S.nparams * log(length(available_case));
cache{j} = add_to_cache(cache{j}, ps, score, max_cache_size, n);
overwrite_count = cache{j}.overwrite_count;
end
function [found, score] = check_cache(cj, ps)
N = size(cj.masks, 2);
mask = false(1, N);
mask(ps) = true;
found = false;
score = 0;
for i = 1:size(cj.masks, 1)
    if isequal(cj.masks(i,:), mask)
        found = true;
        score = cj.scores(i);
        return;
    end
end
end
function cj_new = add_to_cache(cj, ps, score, max_cache_size, total_nodes)
if nargin < 5 || isempty(total_nodes)
    if ~isempty(cj.masks) && size(cj.masks, 2) > 0
        total_nodes = size(cj.masks, 2);
    else
        error('total_nodes is required when the parent-set cache is empty.');
    end
end
mask = false(1, total_nodes);
if ~isempty(ps)
    mask(unique(ps)) = true;
end
if isfield(cj, 'overwrite_count') == false
    cj.overwrite_count = 0;
end
if isempty(cj.masks)
    cj.masks = mask;
    cj.scores = score;
else
    if size(cj.masks, 1) >= max_cache_size
        cj.masks = cj.masks(2:end, :);
        cj.scores = cj.scores(2:end);
        cj.overwrite_count = cj.overwrite_count + 1;
    end
    cj.masks(end+1, :) = mask;
    cj.scores(end+1) = score;
end
cj_new = cj;
end
