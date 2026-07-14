%% AESL-GA HSM analysis for datasets not covered by the default AESL run.
%
% The previous AESL-HSM run covered Alarm and HEPAR. This script completes
% the same analysis on the remaining CS-CGA benchmark networks.

previous_datasets = getenv('AESL_HSM_ANALYSIS_DATASETS');
cleanup_datasets = onCleanup(@() setenv('AESL_HSM_ANALYSIS_DATASETS', previous_datasets));

setenv('AESL_HSM_ANALYSIS_DATASETS', 'Asia,INS,Water,Hailfinder,Win95pts,AND');
Main_AESL_HSM_Analysis;
