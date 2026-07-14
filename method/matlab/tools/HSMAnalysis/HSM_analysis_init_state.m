function state = HSM_analysis_init_state(variant, node_count)
state = struct();
state.variantName = variant.name;
state.mode = variant.mode;
state.nodeCount = node_count;
state.duplicateHistory = containers.Map('KeyType', 'char', 'ValueType', 'logical');
state.globalQueries = 0;
state.globalHits = 0;
state.localQueries = 0;
state.localHits = 0;
state.overwrites = 0;
state.globalKeyChars = 0;
state.localKeyChars = 0;
state.globalCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
state.localHash = containers.Map('KeyType', 'char', 'ValueType', 'any');
state.localCache = [];
state.bntCache = [];
state.parentCacheSize = NaN;
state.bntCacheCapacity = NaN;
state.largeHashBytes = NaN;
switch lower(variant.mode)
    case {'bnt', 'gb'}
        state.bntCacheCapacity = 4 * max(64, node_count) * max(64, node_count);
        state.bntCache = score_init_cache(node_count, state.bntCacheCapacity);
    case 'large-hash'
        if isfield(variant, 'largeHashBytes') && ~isempty(variant.largeHashBytes)
            state.largeHashBytes = variant.largeHashBytes;
        else
            state.largeHashBytes = 4 * 1024^3;
        end
    case 'hsm'
        if isfield(variant, 'parentCacheSize') && ~isempty(variant.parentCacheSize)
            state.parentCacheSize = variant.parentCacheSize;
        elseif isfield(variant, 'parentCacheMultiplier') && ~isempty(variant.parentCacheMultiplier)
            state.parentCacheSize = variant.parentCacheMultiplier * node_count;
        else
            state.parentCacheSize = 200 * node_count;
        end
        state.localCache = cell(1, node_count);
        for node_id = 1:node_count
            state.localCache{node_id} = struct('masks', [], 'scores', [], 'overwrite_count', 0);
        end
    otherwise
        error('HSM_analysis_init_state:UnknownMode', 'Unknown cache mode: %s', variant.mode);
end
end
