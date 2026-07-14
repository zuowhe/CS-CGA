function Bnets = prepare_benchmark_datasets(DsS, trial, Flag, output_dir)
    RootPath = fileparts(mfilename('fullpath'));
    if nargin < 4 || isempty(output_dir)
        output_dir = fullfile(RootPath, 'datasets');
    end
    Bnets = cell(0,0);
    for j = 1:size(DsS,2)
        switch DsS{j}{1,1}
            case 'Asia',        bnet = mk_asia_bnet;
            case 'Cancer',      bnet = mk_cancer_bnet;
            case 'Earthquake',  bnet = mk_earthquake_bnet;
            case 'Sachs',       bnet = mk_sachs_bnet;
            case 'Survey',      bnet = mk_survey_bnet;
            case 'Alarm',       bnet = mk_alarm_bnet;
            case 'Barley',      bnet = mk_barley_bnet;
            case 'INS',         bnet = mk_insur_bnet;
            case 'Mildew',      bnet = mk_mildew_bnet;
            case 'Water',       bnet = mk_water_bnet;
            case 'Hailfinder',  bnet = mk_hailfinder_bnet;
            case 'HEPAR',       bnet = mk_hepar2_bnet;
            case 'Win95pts',    bnet = mk_win95pts_bnet;
            case 'AND',         bnet = mk_andes_bnet;
            case 'Pathfinder',  bnet = mk_pathfinder_bnet;
        end
        Bnets{end+1} = {DsS{j}{1,1},bnet};
        if Flag
            for i = 1:size(DsS{j}{1,2},2)
                    Ds_size = DsS{j}{1,2}(i);
                    str = sprintf('%s%s',DsS{j}{1,1},num2str(Ds_size));
                    fprintf('Generating %s datasets.     Start time [%s]\n',str,datestr(now));
                    [~] = generate_BN_samples(str,Ds_size,trial,bnet,output_dir);
                    fprintf('- Dataset %s is generated. Finish time [%s]\n',str,datestr(now));
            end
        else
            fprintf('No need to generate new datasets.\n')
        end
    end
end
