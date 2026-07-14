function [p] = CS_CGA_random_bitflip_mutation(N,l_map,p)
l_cnt = size(l_map,1);
m = 1/l_cnt;
for l=1:l_cnt
    j = l_map(l,1);     k = l_map(l,2);
    for i=1:N
        l_val = get_allele(p{i}(j,k),p{i}(k,j));
        if m >= rand
            l_val_new = mod(l_val + round(rand),3)+1;
            switch l_val_new
                case 1
                    p{i}(j,k)=false; p{i}(k,j)=false;
                case 2
                    p{i}(j,k)=false; p{i}(k,j)=true;
                case 3
                    p{i}(j,k)=true; p{i}(k,j)=false;
            end
        end
    end
end
end
