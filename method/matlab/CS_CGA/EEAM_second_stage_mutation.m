function [p] = EEAM_second_stage_mutation(l_map,p,all_repeats)
	repeat_size = size(all_repeats,2);
    l_cnt = size(l_map,1);
    nrep = size(all_repeats,2);
    m = 1/l_cnt;
    for L=1:l_cnt
        j = l_map(L,1);     k = l_map(L,2);
        for I=1:nrep
            i = all_repeats(I);
            l_val = get_allele(p{i}(j,k),p{i}(k,j));
            if m >= rand
                switch l_val
                    case 1
                        if rand>0.5
                            p{i}(j,k)=false; p{i}(k,j)=true;
                        else
                            p{i}(j,k)=true; p{i}(k,j)=false;
                        end
                    case 2
                        if l_map(L,3)>rand
                            p{i}(j,k)=true; p{i}(k,j)=false;
                        else
                            p{i}(j,k)=false; p{i}(k,j)=false;
                        end
                    case 3
                        if l_map(L,3)>rand
                            p{i}(j,k)=false; p{i}(k,j)=true;
                        else
                            p{i}(j,k)=false; p{i}(k,j)=false;
                        end
                end
            end
        end
    end
end
