function [MI,norm_MI] = get_MI(data, node_sizes)
bns = size(node_sizes, 2);
m = size(data, 2);
MI = zeros(bns);
for n1 = 1:bns
    for n2 = n1+1 :bns
        cnt1 = zeros(1, node_sizes(n1));
        cnt2 = zeros(1, node_sizes(n2));
        cnt12 = zeros(node_sizes(n1), node_sizes(n2));
        for k = 1:m
            a = data(n1, k);
            b = data(n2, k);
            cnt1(a) = cnt1(a) +1;
            cnt2(b) = cnt2(b) +1;
            cnt12(a, b) = cnt12(a, b) +1;
        end
        cnt1 = cnt1 / m;
        cnt2 = cnt2 / m;
        cnt12 = cnt12 / m;
        for i = 1:node_sizes(n1)
            for j = 1:node_sizes(n2)
                if cnt12(i,j) > 0
                    delta_MI = cnt12(i,j) * log2(cnt12(i,j) / (cnt1(i)*cnt2(j)) );
                    MI(n1, n2) = MI(n1, n2) + delta_MI;
                    MI(n2, n1) = MI(n1, n2);
                end
            end
        end
    end
end
norm_MI = zeros(bns);
MMI = max(MI);
for i = 1:bns-1
    for j = i:bns
        norm_MI(i,j) = max(MI(i,j)/MMI(i), MI(i,j)/MMI(j));
        norm_MI(j,i) = norm_MI(i,j);
    end
end
end
