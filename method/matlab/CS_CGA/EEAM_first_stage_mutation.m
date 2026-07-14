function [p,m1] = EEAM_first_stage_mutation(N,l_map,p,Dif_BIC,M,i)
l_cnt = size(l_map,1);
beta = (M-i)/M;
m1 = calculateMutationRateLog(l_cnt, Dif_BIC, beta);
for l=1:l_cnt
    j = l_map(l,1);     k = l_map(l,2);
    for i=1:N
        l_val = get_allele(p{i}(j,k),p{i}(k,j));
        if m1 >= rand
            switch l_val
                case 1
                    if rand>0.5
                        p{i}(j,k)=false; p{i}(k,j)=true;
                    else
                        p{i}(j,k)=true; p{i}(k,j)=false;
                    end
                case 2
                    if l_map(l,3)>rand
                        p{i}(j,k)=true; p{i}(k,j)=false;
                    else
                        p{i}(j,k)=false; p{i}(k,j)=false;
                    end
                case 3
                    if l_map(l,3)>rand
                        p{i}(j,k)=false; p{i}(k,j)=true;
                    else
                        p{i}(j,k)=false; p{i}(k,j)=false;
                    end
            end
        end
    end
end
end
function mutation_rate = calculateMutationRateLog(n, Dif_BIC, beta)
    % EEAM increases mutation pressure when the current population falls behind the archive best.
    gamma = 1 / (1 + exp(Dif_BIC));
    if isinf(gamma) || isnan(gamma)
        gamma = 0;
    end
    mutation_rate = (1 / n) * (gamma+beta);
end
