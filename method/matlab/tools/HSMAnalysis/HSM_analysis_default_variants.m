function variants = HSM_analysis_default_variants()
variants = struct('name', {}, 'mode', {}, 'parentCacheMultiplier', {}, ...
    'parentCacheSize', {}, 'largeHashBytes', {});
variants(end + 1) = make_variant('BNT-Cache', 'bnt', [], [], []);
variants(end + 1) = make_variant('GB-Cache', 'gb', [], [], []);
variants(end + 1) = make_variant('Large-Hash', 'large-hash', [], [], 4 * 1024^3);
variants(end + 1) = make_variant('HSM-50n', 'hsm', 50, [], []);
variants(end + 1) = make_variant('HSM-100n', 'hsm', 100, [], []);
variants(end + 1) = make_variant('HSM-200n', 'hsm', 200, [], []);
variants(end + 1) = make_variant('HSM-500n', 'hsm', 500, [], []);
variants(end + 1) = make_variant('HSM-large', 'hsm', [], 1e9, []);
end
function variant = make_variant(name, mode, parent_cache_multiplier, parent_cache_size, large_hash_bytes)
variant = struct();
variant.name = name;
variant.mode = mode;
variant.parentCacheMultiplier = parent_cache_multiplier;
variant.parentCacheSize = parent_cache_size;
variant.largeHashBytes = large_hash_bytes;
end
