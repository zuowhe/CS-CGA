dbstop if error
setup_CS_CGA_paths;
data_size = [500, 1000, 3000];
DsS = cell(0, 0);
DsS{end+1} = {'Asia', data_size};
DsS{end+1} = {'INS', data_size};
DsS{end+1} = {'Water', data_size};
DsS{end+1} = {'Alarm', data_size};
DsS{end+1} = {'Hailfinder', data_size};
DsS{end+1} = {'HEPAR', data_size};
DsS{end+1} = {'Win95pts', data_size};
DsS{end+1} = {'AND', data_size};
Train_Num = 10;
FlagNewdata = false;
Bnets = prepare_benchmark_datasets(DsS, Train_Num, FlagNewdata);
N = 100;
M = 200;
MP = 7;
tour = 2;
scoring_fn = 'bic';
AblationVariants = {
    'SR-Fix', ...
    'SR-PCR', ...
    'SM-Fix', ...
    'SM-PCR', ...
    'DR-PCR', ...
    'CS-CGA'
};
for a = 1:numel(AblationVariants)
    variant = AblationVariants{a};
    for j = 1:size(DsS, 2)
        BN_Name = DsS{j}{1, 1};
        Ds_set = DsS{j}{1, 2};
        bnet = Bnets{j}{1, 2};
        for i = 1:size(Ds_set, 2)
            Ds = Ds_set(i);
            run_CS_CGA_ablation_on_dataset(variant, Ds, BN_Name, N, M, MP, tour, bnet, Train_Num, scoring_fn);
        end
    end
end
