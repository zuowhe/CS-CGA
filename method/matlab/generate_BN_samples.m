function [data] = generate_BN_samples(str,D,T,bnet,RootPath)
n = size(bnet.dag,1);
    data = cell(1,T);
for t=1:T
    data{t} = zeros(n,D);
    for i=1:D
        sample = sample_bnet(bnet);
        data{t}(:,i) = [sample{:}];
    end
end
eval([str, '=data;']);
save(fullfile(RootPath,'datasets', [str '.mat']), str);
end
