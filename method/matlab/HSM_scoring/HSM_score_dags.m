function [score, cache, total_overwrite] = HSM_score_dags(data, ns, dags, varargin)
[n, ncases] = size(data);
type = cell(1, n);
params = cell(1, n);
for i = 1:n
    type{i} = 'tabular';
    params{i} = {'prior_type', 'dirichlet', 'dirichlet_weight', 1};
end
scoring_fn = 'bayesian';
discrete = 1:n;
clamped = zeros(n, ncases);
cache = [];
max_cache_size = 100;
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}
        case 'scoring_fn', scoring_fn = args{i+1};
        case 'type',       type = args{i+1};
        case 'discrete',   discrete = args{i+1};
        case 'clamped',    clamped = args{i+1};
        case 'params',     params = args{i+1};
        case 'cache',      cache = args{i+1};
        case 'max_cache_size', max_cache_size = args{i+1};
    end
end
NG = length(dags);
score = zeros(1, NG);
total_overwrite = 0;
for g = 1:NG
    if isempty(dags{g})
        score(g) = -Inf;
        continue;
    end
    for j = 1:n
        ps = parents(dags{g}, j);
        [scor, cache, hit_count, overwrite_count] = HSM_score_family(j, ps, type{j}, ns, discrete, data, params{j}, cache, max_cache_size);
        score(g) = score(g) + scor;
        total_overwrite = total_overwrite + overwrite_count;
    end
end
end
