function stats = HSM_analysis_cache_stats(state)
stats = struct();
stats.variant = state.variantName;
stats.cacheMode = state.mode;
stats.globalQueries = state.globalQueries;
stats.globalHits = state.globalHits;
stats.localQueries = state.localQueries;
stats.localHits = state.localHits;
stats.overwrites = state.overwrites;
stats.globalEntries = 0;
stats.localEntries = 0;
stats.bntCacheCapacity = state.bntCacheCapacity;
stats.parentCacheSize = state.parentCacheSize;
stats.largeHashTargetMB = bytes_to_mb(state.largeHashBytes);
stats.estimatedCacheMB = 0;
stats.matlabUsedMB = get_matlab_used_mb();
if state.globalQueries > 0
    stats.globalHitRate = state.globalHits / state.globalQueries;
else
    stats.globalHitRate = NaN;
end
if state.localQueries > 0
    stats.localHitRate = state.localHits / state.localQueries;
else
    stats.localHitRate = NaN;
end
switch lower(state.mode)
    case 'bnt'
        stats.localEntries = count_bnt_entries(state.bntCache);
        stats.estimatedCacheMB = bytes_to_mb(numel(state.bntCache) * 8);
    case 'gb'
        stats.globalEntries = state.globalCache.Count;
        stats.localEntries = count_bnt_entries(state.bntCache);
        bnt_bytes = numel(state.bntCache) * 8;
        global_bytes = state.globalKeyChars + 8 * state.globalCache.Count;
        stats.estimatedCacheMB = bytes_to_mb(bnt_bytes + global_bytes);
    case 'large-hash'
        stats.globalEntries = state.globalCache.Count;
        stats.localEntries = state.localHash.Count;
        global_bytes = state.globalKeyChars + 8 * state.globalCache.Count;
        local_bytes = state.localKeyChars + 8 * state.localHash.Count;
        stats.estimatedCacheMB = bytes_to_mb(global_bytes + local_bytes);
    case 'hsm'
        stats.globalEntries = state.globalCache.Count;
        stats.localEntries = count_hsm_entries(state.localCache);
        global_bytes = state.globalKeyChars + 8 * state.globalCache.Count;
        local_bytes = stats.localEntries * (state.nodeCount + 8);
        stats.estimatedCacheMB = bytes_to_mb(global_bytes + local_bytes);
    otherwise
        error('HSM_analysis_cache_stats:UnknownMode', 'Unknown cache mode: %s', state.mode);
end
end
function entries = count_bnt_entries(cache)
entries = 0;
if isempty(cache)
    return
end
node_count = size(cache, 2) - 3;
if node_count <= 0 || size(cache, 1) < 2
    return
end
entries = sum(cache(2:end, node_count + 1) > 0);
end
function entries = count_hsm_entries(local_cache)
entries = 0;
if isempty(local_cache)
    return
end
for node_id = 1:numel(local_cache)
    if isfield(local_cache{node_id}, 'masks') && ~isempty(local_cache{node_id}.masks)
        entries = entries + size(local_cache{node_id}.masks, 1);
    end
end
end
function mb = bytes_to_mb(bytes)
if isnan(bytes)
    mb = NaN;
else
    mb = bytes / 1024^2;
end
end
function used_mb = get_matlab_used_mb()
used_mb = NaN;
try
    mem_info = memory;
    if isfield(mem_info, 'MemUsedMATLAB')
        used_mb = mem_info.MemUsedMATLAB / 1024^2;
    end
catch
end
end
