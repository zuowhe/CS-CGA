function [] = dag2file(dag, path, algo, net, Ds, N, M, T)
saves_filename_dag = sprintf('%s/DAG_%s%s_%s_%s_%s-%s.csv',...
    path,net,num2str(Ds),num2str(N),num2str(M),algo,num2str(T,'%02d'));
save_file = fopen(saves_filename_dag, 'a+');
data = dag{1};
bns = size(data,1);
for i = 1:bns
    for j = 1:bns
        fprintf(save_file,'%1d,',data(i,j));
    end
    fprintf(save_file,'\n');
end
fprintf(save_file,'\n');
fclose(save_file);
end
