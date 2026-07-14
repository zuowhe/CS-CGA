function [conf, conf_norm] = get_confidence(pop)
N = max(size(pop,2),size(pop,1));
bns = size(pop{1},1);
conf = zeros(bns);
for i = 1:bns
    for j = i:bns
        for ii = 1:N
            conf(i,j) = conf(i,j) + pop{ii}(i,j);
            conf(j,i) = conf(j,i) + pop{ii}(j,i);
        end
    end
end
[~,~,conf_list] = find(conf);
eps = 0.0000001;
conf_min = min(conf_list);
conf_max = max(conf_list);
conf_range = conf_max - conf_min;
conf_norm = (conf-conf_min)/(eps+conf_range);
conf_norm(conf_norm < 0) = 0;
end
