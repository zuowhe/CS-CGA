function setup_HSM_analysis_paths(root_dir)
if nargin < 1 || isempty(root_dir)
    root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
release_root = fileparts(fileparts(root_dir));
bnt_dir = fullfile(release_root, 'third_party', 'bnt-master');
addpath(root_dir, '-begin');
addpath(fullfile(root_dir, 'tools', 'HSMAnalysis'), '-begin');
addpath(fullfile(root_dir, 'CS_CGA'), '-begin');
addpath(fullfile(root_dir, 'HSM_scoring'), '-begin');
addpath(fullfile(root_dir, 'Subprogram'), '-begin');
addpath(fullfile(root_dir, 'datasets'), '-begin');
addpath(genpath(fullfile(root_dir, 'Subprogram', 'Bif2Bnt')));
if exist(bnt_dir, 'dir')
    addpath(genpath(bnt_dir));
else
    addpath(genpath(fullfile(root_dir, 'bnt-master')));
end
end
