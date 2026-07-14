function cfg = CS_CGA_ablation_config(variant)
key = upper(strrep(strrep(strtrim(variant), '_', '-'), ' ', '-'));
cfg = struct();
cfg.key = key;
cfg.fixed_alpha = 0.01;
cfg.use_pcr = false;
cfg.use_mi_mutation = false;
cfg.use_second_stage = false;
cfg.use_second_stage_mi = false;
cfg.output_name = '';
cfg.description = '';
switch key
    case 'SR-FIX'
        cfg.output_name = 'SR-Fix';
        cfg.description = 'Single-stage random mutation with fixed alpha=0.01.';
    case 'SR-PCR'
        cfg.output_name = 'SR-PCR';
        cfg.use_pcr = true;
        cfg.description = 'Single-stage random mutation with progressive constraint relaxation.';
    case 'SM-FIX'
        cfg.output_name = 'SM-Fix';
        cfg.use_mi_mutation = true;
        cfg.description = 'Single-stage MI-guided mutation with fixed alpha=0.01.';
    case 'SM-PCR'
        cfg.output_name = 'SM-PCR';
        cfg.use_pcr = true;
        cfg.use_mi_mutation = true;
        cfg.description = 'Single-stage MI-guided mutation with progressive constraint relaxation.';
    case 'DR-PCR'
        cfg.output_name = 'DR-PCR';
        cfg.use_pcr = true;
        cfg.use_second_stage = true;
        cfg.description = 'Dual-stage random mutation with progressive constraint relaxation.';
    case {'CS-CGA', 'CS_CGA'}
        cfg.key = 'CS-CGA';
        cfg.output_name = 'CS-CGA';
        cfg.use_pcr = true;
        cfg.use_mi_mutation = true;
        cfg.use_second_stage = true;
        cfg.use_second_stage_mi = true;
        cfg.description = 'Full method with PCR, EEAM, MI guidance, and HSM scoring.';
    otherwise
        error('Unknown ablation variant: %s', variant);
end
end
