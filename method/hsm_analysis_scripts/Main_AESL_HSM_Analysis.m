%% AESL-GA cache behavior and cache-size tradeoff analysis.
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'tools', 'AESLHSMAnalysis'), '-begin');
setup_AESL_HSM_analysis_paths(root_dir);

dataset_override = strtrim(getenv('AESL_HSM_ANALYSIS_DATASETS'));
if isempty(dataset_override)
    datasets = {'Alarm', 'HEPAR'};
else
    datasets = strtrim(strsplit(dataset_override, ','));
end

sample_size_override = strtrim(getenv('AESL_HSM_ANALYSIS_SAMPLE_SIZE'));
if isempty(sample_size_override)
    sample_size = 1000;
else
    sample_size = str2double(sample_size_override);
end

train_index = 1;
train_num = 1;
flag_new_data = false;

N = 100;
M = 200;
MP = 7;
scoring_fn = 'bic';
tol = 0.01;
alpha = 0.9;
d = struct('m0', 1 / 5, 'M0', 3 / 5, 'm1', 1 / 10, 'M1', 1 / 2);
base_seed = 20260627;
timeout_override = strtrim(getenv('AESL_HSM_ANALYSIS_TIMEOUT_SEC'));
if isempty(timeout_override)
    time_limit_sec = 10000;
else
    time_limit_sec = str2double(timeout_override);
end

DsS = cell(0, 0);
for dataset_id = 1:numel(datasets)
    DsS{end + 1} = {datasets{dataset_id}, sample_size}; %#ok<SAGROW>
end
Bnets = Generate_dataset(DsS, train_num, flag_new_data);

variants = HSM_analysis_default_variants();
variant_override = strtrim(getenv('AESL_HSM_ANALYSIS_VARIANTS'));
if ~isempty(variant_override)
    requested_variants = strtrim(strsplit(variant_override, ','));
    available_names = {variants.name};
    missing_variants = requested_variants(~ismember(requested_variants, available_names));
    if ~isempty(missing_variants)
        error('Main_AESL_HSM_Analysis:UnknownVariant', ...
            'Unknown cache variant(s): %s', strjoin(missing_variants, ', '));
    end
    variants = variants(ismember(available_names, requested_variants));
end

run_id = datestr(now, 'yyyymmdd_HHMMSS');
result_dir = fullfile(root_dir, 'results', ['AESL_HSM_analysis_', run_id]);
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

summary_file = fullfile(result_dir, 'aesl_hsm_analysis_summary.csv');
write_aesl_hsm_analysis_header(summary_file);

fprintf('AESL-HSM analysis output: %s\n', summary_file);

for dataset_id = 1:numel(datasets)
    BN_Name = datasets{dataset_id};
    bnet = Bnets{dataset_id}{1, 2};
    bnet.dag = logical(bnet.dag);
    data_name = sprintf('%s%d', BN_Name, sample_size);
    loaded = load(data_name);
    TrainData = loaded.(data_name);
    data = TrainData{train_index};
    super_structure = get_CI_test(data, bnet, tol);

    for variant_id = 1:numel(variants)
        variant = variants(variant_id);
        rng(base_seed + dataset_id, 'twister');
        fprintf('[%s] Dataset=%s-%d Algorithm=AESL-GA Variant=%s\n', ...
            datestr(now), BN_Name, sample_size, variant.name);
        t_start = tic;

        [dag, best_score, ~, iterations, cache_stats] = AESL_HSM_analysis_process( ...
            super_structure, data, N, M, MP, alpha, d, scoring_fn, bnet, variant, time_limit_sec);

        runtime_sec = toc(t_start);
        [f1, se, sp, precision, shd, TP, TP2, FN, FP, TN] = eval_dags_adjust(dag, bnet.dag, 1);

        row = struct();
        row.algorithm = 'AESL-GA';
        row.dataset = BN_Name;
        row.sampleSize = sample_size;
        row.variant = variant.name;
        row.status = cache_stats.status;
        row.timeLimitSec = cache_stats.timeLimitSec;
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
        row.localEntries = cache_stats.localEntries;
        row.overwrites = cache_stats.overwrites;
        row.bntCacheCapacity = cache_stats.bntCacheCapacity;
        row.parentCacheSize = cache_stats.parentCacheSize;
        row.largeHashTargetMB = cache_stats.largeHashTargetMB;
        row.estimatedCacheKB = cache_stats.estimatedCacheKB;
        row.estimatedCacheMB = cache_stats.estimatedCacheMB;
        row.globalCacheKB = cache_stats.globalCacheKB;
        row.localCacheKB = cache_stats.localCacheKB;
        row.totalCacheKB = cache_stats.totalCacheKB;
        row.globalCacheMB = cache_stats.globalCacheMB;
        row.localCacheMB = cache_stats.localCacheMB;
        row.totalCacheMB = cache_stats.totalCacheMB;
        row.localLimitEntries = cache_stats.localLimitEntries;
        row.finalLocalEntries = cache_stats.finalLocalEntries;
        row.peakLocalEntries = cache_stats.peakLocalEntries;
        row.localLimitKB = cache_stats.localLimitKB;
        row.localUsageRatio = cache_stats.localUsageRatio;
        row.matlabUsedMB = cache_stats.matlabUsedMB;
        row.initialAlpha = cache_stats.initialAlpha;
        row.finalAlpha = cache_stats.finalAlpha;
        row.alphaUpdates = cache_stats.alphaUpdates;
        row.searchSpaceEdges = cache_stats.searchSpaceEdges;
        row.TP = TP;
        row.TP2 = TP2;
        row.FN = FN;
        row.FP = FP;
        row.TN = TN;

        print_aesl_hsm_analysis_row(row);
        append_aesl_hsm_analysis_row(summary_file, row);
    end
end

fprintf('AESL-HSM analysis finished: %s\n', summary_file);

function write_aesl_hsm_analysis_header(filename)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['algorithm,dataset,sampleSize,variant,cacheMode,nodes,iterations,f1,sensitivity,' ...
    'specificity,precision,shd,bestScore,runtimeSec,scoringTimeSec,scoringFraction,' ...
    'globalQueries,globalHits,' ...
    'globalHitRate,localQueries,localHits,localHitRate,globalEntries,localEntries,' ...
    'overwrites,bntCacheCapacity,parentCacheSize,largeHashTargetMB,estimatedCacheKB,' ...
    'estimatedCacheMB,globalCacheKB,localCacheKB,totalCacheKB,globalCacheMB,localCacheMB,' ...
    'totalCacheMB,localLimitEntries,finalLocalEntries,peakLocalEntries,localLimitKB,localUsageRatio,' ...
    'matlabUsedMB,initialAlpha,finalAlpha,alphaUpdates,searchSpaceEdges,TP,TP2,FN,FP,TN,status,timeLimitSec\n']);
end

function append_aesl_hsm_analysis_row(filename, row)
fid = fopen(filename, 'a');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, ['%s,%s,%d,%s,%s,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,' ...
    '%d,%d,%.8g,%d,%d,%.8g,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,' ...
    '%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,' ...
    '%.8g,%.8g,%d,%d,%d,%d,%d,%d,%d,%s,%.8g\n'], ...
    row.algorithm, row.dataset, row.sampleSize, row.variant, row.cacheMode, row.nodes, row.iterations, ...
    row.f1, row.sensitivity, row.specificity, row.precision, row.shd, row.bestScore, ...
    row.runtimeSec, row.scoringTimeSec, row.scoringFraction, ...
    row.globalQueries, row.globalHits, row.globalHitRate, ...
    row.localQueries, row.localHits, row.localHitRate, row.globalEntries, row.localEntries, ...
    row.overwrites, row.bntCacheCapacity, row.parentCacheSize, row.largeHashTargetMB, ...
    row.estimatedCacheKB, row.estimatedCacheMB, row.globalCacheKB, row.localCacheKB, ...
    row.totalCacheKB, row.globalCacheMB, row.localCacheMB, row.totalCacheMB, ...
    row.localLimitEntries, row.finalLocalEntries, row.peakLocalEntries, row.localLimitKB, ...
    row.localUsageRatio, row.matlabUsedMB, ...
    row.initialAlpha, row.finalAlpha, ...
    row.alphaUpdates, row.searchSpaceEdges, row.TP, row.TP2, row.FN, row.FP, row.TN, ...
    row.status, row.timeLimitSec);
end

function print_aesl_hsm_analysis_row(row)
fprintf(['RESULT dataset=%s variant=%s runtime=%.3fs globalHit=%.5f ' ...
    'localHit=%.5f overwrites=%d globalKB=%.3f localKB=%.3f totalKB=%.3f status=%s\n'], ...
    row.dataset, row.variant, row.runtimeSec, row.globalHitRate, ...
    row.localHitRate, row.overwrites, row.globalCacheKB, row.localCacheKB, row.totalCacheKB, row.status);
end
