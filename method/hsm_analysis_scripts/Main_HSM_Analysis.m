%% HSM cache behavior and cache-size tradeoff analysis.
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'tools', 'HSMAnalysis'), '-begin');
setup_HSM_analysis_paths(root_dir);

dataset_override = strtrim(getenv('HSM_ANALYSIS_DATASETS'));
if isempty(dataset_override)
    datasets = {'Asia', 'INS', 'Water', 'Alarm', 'Hailfinder', 'HEPAR', 'Win95pts', 'AND'};
else
    datasets = strtrim(strsplit(dataset_override, ','));
end
sample_size = 1000;
train_index = 1;
train_num = 1;
flag_new_data = false;

N = 100;
M = 200;
MP = 7;
tour = 2;
scoring_fn = 'bic';
ci_update_rule = 'paper';
base_seed = 20260609;

DsS = cell(0, 0);
for dataset_id = 1:numel(datasets)
    DsS{end + 1} = {datasets{dataset_id}, sample_size}; %#ok<SAGROW>
end
Bnets = Generate_dataset(DsS, train_num, flag_new_data);

variants = HSM_analysis_default_variants();
variants = add_cscga_only_variants(variants);
variant_override = strtrim(getenv('HSM_ANALYSIS_VARIANTS'));
if ~isempty(variant_override)
    requested_variants = strtrim(strsplit(variant_override, ','));
    available_names = {variants.name};
    missing_variants = requested_variants(~ismember(requested_variants, available_names));
    if ~isempty(missing_variants)
        error('Main_HSM_Analysis:UnknownVariant', ...
            'Unknown cache variant(s): %s', strjoin(missing_variants, ', '));
    end
    variants = variants(ismember(available_names, requested_variants));
end

run_id = datestr(now, 'yyyymmdd_HHMMSS');
result_dir = fullfile(root_dir, 'results', ['HSM_analysis_', run_id]);
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

summary_file = fullfile(result_dir, 'hsm_analysis_summary.csv');
write_hsm_analysis_header(summary_file);

fprintf('HSM analysis output: %s\n', summary_file);

for dataset_id = 1:numel(datasets)
    BN_Name = datasets{dataset_id};
    bnet = Bnets{dataset_id}{1, 2};
    data_name = sprintf('%s%d', BN_Name, sample_size);
    loaded = load(data_name);
    TrainData = loaded.(data_name);
    data = TrainData{train_index};
    p_value = Pvalue_get_CI_test(data, bnet, 0.01);

    for variant_id = 1:numel(variants)
        variant = variants(variant_id);
        rng(base_seed + dataset_id, 'twister');
        fprintf('[%s] Dataset=%s-%d Variant=%s\n', datestr(now), BN_Name, sample_size, variant.name);
        t_start = tic;

        [dag, best_score, ~, iterations, cache_stats] = HSM_analysis_process( ...
            data, N, M, MP, scoring_fn, bnet, tour, p_value, variant, ci_update_rule);

        runtime_sec = toc(t_start);
        [f1, se, sp, precision, shd, TP, TP2, FN, FP, TN] = eval_dags_adjust(dag, bnet.dag, 1);

        row = struct();
        row.dataset = BN_Name;
        row.sampleSize = sample_size;
        row.variant = variant.name;
        row.cacheMode = cache_stats.cacheMode;
        row.nodes = size(bnet.dag, 1);
        row.iterations = iterations;
        row.f1 = f1;
        row.sensitivity = se;
        row.specificity = sp;
        row.precision = precision;
        row.shd = shd;
        row.bestScore = best_score;
        row.runtimeSec = runtime_sec;
        row.scoringTimeSec = cache_stats.scoringTimeSec;
        row.scoringFraction = cache_stats.scoringTimeSec / runtime_sec;
        row.globalQueries = cache_stats.globalQueries;
        row.globalHits = cache_stats.globalHits;
        row.globalHitRate = cache_stats.globalHitRate;
        row.localQueries = cache_stats.localQueries;
        row.localHits = cache_stats.localHits;
        row.localHitRate = cache_stats.localHitRate;
        row.globalEntries = cache_stats.globalEntries;
        row.parentSetCacheEntries = cache_stats.finalLocalEntries;
        row.overwrites = cache_stats.overwrites;
        row.parentCacheSize = cache_stats.parentCacheSize;
        row.localLimitEntries = cache_stats.localLimitEntries;
        row.peakLocalEntries = cache_stats.peakLocalEntries;
        row.localUsageRatio = cache_stats.localUsageRatio;
        row.initialAlpha = cache_stats.initialAlpha;
        row.finalAlpha = cache_stats.finalAlpha;
        row.alphaUpdates = cache_stats.alphaUpdates;
        row.searchSpaceEdges = cache_stats.searchSpaceEdges;
        row.TP = TP;
        row.TP2 = TP2;
        row.FN = FN;
        row.FP = FP;
        row.TN = TN;

        print_hsm_analysis_row(row);
        append_hsm_analysis_row(summary_file, row);
    end
end

fprintf('HSM analysis finished: %s\n', summary_file);

function write_hsm_analysis_header(filename)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['dataset,sampleSize,variant,cacheMode,nodes,iterations,f1,sensitivity,' ...
    'specificity,precision,shd,bestScore,runtimeSec,scoringTimeSec,scoringFraction,' ...
    'globalQueries,globalHits,' ...
    'globalHitRate,localQueries,localHits,localHitRate,globalEntries,parentSetCacheEntries,' ...
    'overwrites,parentCacheSize,localLimitEntries,peakLocalEntries,localUsageRatio,' ...
    'initialAlpha,finalAlpha,alphaUpdates,searchSpaceEdges,TP,TP2,FN,FP,TN\n']);
end

function append_hsm_analysis_row(filename, row)
fid = fopen(filename, 'a');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['%s,%d,%s,%s,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,' ...
    '%d,%d,%.8g,%d,%d,%.8g,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,' ...
    '%.8g,%.8g,%d,%d,%d,%d,%d,%d,%d\n'], ...
    row.dataset, row.sampleSize, row.variant, row.cacheMode, row.nodes, row.iterations, ...
    row.f1, row.sensitivity, row.specificity, row.precision, row.shd, row.bestScore, ...
    row.runtimeSec, row.scoringTimeSec, row.scoringFraction, ...
    row.globalQueries, row.globalHits, row.globalHitRate, ...
    row.localQueries, row.localHits, row.localHitRate, row.globalEntries, ...
    row.parentSetCacheEntries, row.overwrites, row.parentCacheSize, ...
    row.localLimitEntries, row.peakLocalEntries, row.localUsageRatio, ...
    row.initialAlpha, row.finalAlpha, row.alphaUpdates, row.searchSpaceEdges, ...
    row.TP, row.TP2, row.FN, row.FP, row.TN);
end

function print_hsm_analysis_row(row)
fprintf(['RESULT dataset=%s variant=%s runtime=%.3fs globalHit=%.5f ' ...
    'localHit=%.5f overwrites=%d parentSetEntries=%d limitEntries=%.8g\n'], ...
    row.dataset, row.variant, row.runtimeSec, row.globalHitRate, ...
    row.localHitRate, row.overwrites, row.parentSetCacheEntries, row.localLimitEntries);
end

function variants = add_cscga_only_variants(variants)
parent_hash = struct();
parent_hash.name = 'ParentHash-paper';
parent_hash.mode = 'large-hash';
parent_hash.parentCacheMultiplier = [];
parent_hash.parentCacheSize = [];
parent_hash.largeHashBytes = [];
variants(end + 1) = parent_hash;
end
