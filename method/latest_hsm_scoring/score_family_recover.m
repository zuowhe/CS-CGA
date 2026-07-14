function [score, cache, cache_hit_count, overwrite_count] = score_family_recover(j, ps, node_type, ns, discrete, data, args, cache, max_cache_size)
% SCORE_FAMILY_RECOVER 计算节点 j 的 BIC 分数，并使用带容量限制的缓存机制
% 输出：
%   score           - 当前父集对应的评分
%   cache           - 更新后的缓存结构
%   cache_hit_count - 本次调用中缓存命中的次数（0 或 1）
%   overwrite_count - 在本次调用中发生的覆盖次数

if nargin < 9 || isempty(max_cache_size)
    max_cache_size = 100;
end

if nargin < 9 || isempty(cache)
    n_nodes = size(data, 1);
    cache = cell(1, n_nodes);
    for k = 1:n_nodes
        cache{k} = struct('masks', [], 'scores', []);
    end
end

% 初始化命中和覆盖计数器
cache_hit_count = 0;
overwrite_count = 0;

misv = -9999;
ccc = iscell(data);
if ccc
    data = bnt_to_mat(data, misv);
end

ps = unique(ps);
[n, ncases] = size(data);

% 查找缓存
[found, score] = check_cache(cache{j}, ps);
if found
    cache_hit_count = 1; % 命中一次
    return;
end

% 构造子图
dag = zeros(n);
if ~isempty(ps)
    dag(ps, j) = 1;
    ps = sort(ps);
end

% 创建网络和 CPD
bnet = mk_bnet(dag, ns, 'discrete', discrete);
fname = sprintf('%s_CPD', node_type);
if isempty(args)
    bnet.CPD{j} = feval(fname, bnet, j);
else
    bnet.CPD{j} = feval(fname, bnet, j, args{:});
end

% 学习参数 & 评分
fam = [ps j];
[~, available_case] = find(data(fam,:) == misv);
available_case = setdiff(1:ncases, available_case);

bnet.CPD{j} = learn_params(bnet.CPD{j}, fam, data(:, available_case), ns, bnet.cnodes);
L = log_prob_node(bnet.CPD{j}, data(j, available_case), data(ps, available_case));
S = struct(bnet.CPD{j});
score = L - 0.5 * S.nparams * log(length(available_case));

% 添加到缓存
if isfield(cache{j}, 'overwrite_count')
    before_overwrite_count = cache{j}.overwrite_count;
else
    before_overwrite_count = 0;
end
cache{j} = add_to_cache(cache{j}, ps, score, max_cache_size, n);

% 获取本次调用新增的覆盖次数
overwrite_count = cache{j}.overwrite_count - before_overwrite_count;

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
% ADD_TO_CACHE 将新的父集及其得分加入缓存，并限制最大条目数，支持循环覆盖
% 
% 输入参数：
%   cj           - 当前缓存结构体
%   ps           - 父节点索引列表（可能为空）
%   score        - 分数
%   max_cache_size - 最大缓存条目数
%   total_nodes  - 总节点数量（用于确定 mask 维度）

% 如果没有传入 total_nodes，则尝试从当前 masks 推断
if nargin < 5 || isempty(total_nodes)
    if ~isempty(cj.masks) && size(cj.masks, 2) > 0
        total_nodes = size(cj.masks, 2);
    else
        error('必须提供 total_nodes 参数或已有缓存数据');
    end
end

% 构造 mask：total_nodes 维度的逻辑向量
mask = false(1, total_nodes);
if ~isempty(ps)
    mask(unique(ps)) = true;
end

% 初始化覆盖计数器（如果不存在）
if isfield(cj, 'overwrite_count') == false
    cj.overwrite_count = 0;
end

% 如果缓存为空，直接初始化
if isempty(cj.masks)
    cj.masks = mask;
    cj.scores = score;
else
    % 超出容量则移除最早记录，并增加覆盖计数
    if size(cj.masks, 1) >= max_cache_size
        cj.masks = cj.masks(2:end, :);
        cj.scores = cj.scores(2:end);
        cj.overwrite_count = cj.overwrite_count + 1; % 增加覆盖计数
    end
    
    % 追加新记录
    cj.masks(end+1, :) = mask;
    cj.scores(end+1) = score;
end

cj_new = cj;
end
