function run_CS_CGA_on_dataset(Ds,BN_Name,N,M,MP,tour,bnet,trial,scoring_fn,ci_update_rule)
AlgoName = 'CS-CGA';
if nargin < 10 || isempty(ci_update_rule)
    ci_update_rule = 'paper';
end
fprintf('Running... %s_%s%s - N=%s M=%s [%s]\n',...
    AlgoName,BN_Name,num2str(Ds),num2str(N),num2str(M),datestr(now));
fprintf('Iter   F1 Score   Sensitivity   Specificity   Precision   HD    BIC Score      Runtime  #Generations TP  TP2  FN  FP  TN\n');
tol = 0.01;
eoer = cell(trial,1);
conv(1:trial) = struct('f1',zeros(1,M),'se',zeros(1,M),'sp',zeros(1,M),'sc',zeros(1,M));
NetSizeStr = sprintf('%s%s',BN_Name,num2str(Ds));
data = load(NetSizeStr);
TrainData = data.(NetSizeStr);
format shortG
NowDate = [datestr(now,10),datestr(now,5),datestr(now,7)];
result_dir = ['results/',NowDate];
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end
saves_filename_result = sprintf('results/%s/%s%s_%s_%s_%s.csv',...
    NowDate,BN_Name,num2str(Ds),num2str(N),num2str(M),AlgoName);
convergence_dir = ['results/',NowDate,'/Convergence_Behavior'];
if ~exist(convergence_dir,'dir')
    mkdir(convergence_dir);
end
for T = 1:trial
    saves_result_file = sprintf('results/%s/Convergence_Behavior/%s%s_%s_%s_%s-%s.csv',...
        NowDate,BN_Name,num2str(Ds),num2str(N),num2str(M),AlgoName,num2str(T,'%02d'));
    eoer{T} = -1*ones(1,9);
    tStart = tic;
    p_value = PCR_compute_pvalue_matrix(TrainData{T},bnet,tol);
    [dag,eoer{T}(1,4),conv(T),eoer{T}(1,6)] = ...
        CS_CGA_search(TrainData{T},N,M,MP,scoring_fn,bnet,tour,saves_result_file,p_value,ci_update_rule);
    eoer{T}(1,5) = toc(tStart);
    [eoer{T}(1,1),eoer{T}(1,2),eoer{T}(1,3),eoer{T}(1,8),eoer{T}(1,7),TP,TP2,FN,FP,TN]= eval_dags_adjust(dag,bnet.dag,1);
    fprintf('%4d  %9.5f    %9.5f    %9.5f  %9.5f  %3.1f   %11.3f  %11.3f      %4d       %4d%4d%4d%4d%4d\n',...
        T,eoer{T}(1),eoer{T}(2),eoer{T}(3),eoer{T}(8),eoer{T}(7),eoer{T}(4),eoer{T}(5),eoer{T}(6),TP,TP2,FN,FP,TN);
    saved_file = fopen(saves_filename_result,'a+');
    fprintf(saved_file,'%4d,%9.5f,%9.5f,%9.5f,%9.5f,%3.1f,%11.3f,%11.3f,%4d,%4d,%4d,%4d,%4d,%4d\n',...
        T,eoer{T}(1),eoer{T}(2),eoer{T}(3),eoer{T}(8),eoer{T}(7),eoer{T}(4),eoer{T}(5),eoer{T}(6),TP,TP2,FN,FP,TN);
    fclose(saved_file);
end
[eoer_avg, eoer_std] = avg_std(eoer,trial);
print_avg_std(eoer_avg,eoer_std);
saved_file = fopen(saves_filename_result,'a+');
fprintf(saved_file,'Avg,%9.5f,%9.5f,%9.5f,%9.5f,%5.3f,%11.3f,%11.3f,%7.3f,%11.3f\n',...
        eoer_avg(1),eoer_avg(2),eoer_avg(3),eoer_avg(8),eoer_avg(7),eoer_avg(4),eoer_avg(5),eoer_avg(6),eoer_avg(9));
fprintf(saved_file,'Std,%9.5f,%9.5f,%9.5f,%9.5f,%5.3f,%11.3f,%11.3f,%7.3f,%11.3f\n',...
        eoer_std(1),eoer_std(2),eoer_std(3),eoer_std(8),eoer_std(7),eoer_std(4),eoer_std(5),eoer_std(6),eoer_std(9));
fclose(saved_file);
end
