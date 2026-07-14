function setup_CS_CGA_paths()
root_dir = fileparts(mfilename('fullpath'));
release_root = fileparts(fileparts(root_dir));
bnt_dir = fullfile(release_root, 'third_party', 'bnt-master');
addpath(root_dir, '-begin');
addpath(fullfile(root_dir, 'CS_CGA'), '-begin');
addpath(fullfile(root_dir, 'HSM_scoring'), '-begin');
addpath(fullfile(root_dir, 'Subprogram'), '-begin');
addpath(genpath(fullfile(root_dir, 'Subprogram', 'Bif2Bnt')));
addpath(fullfile(root_dir, 'tools', 'HSMAnalysis'));
if exist(bnt_dir, 'dir')
    addpath(genpath(bnt_dir));
end
end
