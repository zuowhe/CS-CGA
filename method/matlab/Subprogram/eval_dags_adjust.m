function [f1,se,sp,pr,shd,TP,TP2,FN,FP,TN] = eval_dags_adjust(dags,dag0,N)
bns = size(dag0,1);
totP0 = 0;
totN0 = 0;
for j = 1:bns
    for k = j:bns
        if xor(dag0(j,k),dag0(k,j))
            totP0 = totP0+1;
        else
            totN0 = totN0+1;
        end
    end
end
f1 = zeros(1,N);
se = zeros(1,N);
sp = zeros(1,N);
pr = zeros(1,N);
shd = zeros(1,N);
for i = 1:N
    TP = 0;
    TP2 = 0;
    FP = 0;
    TN = 0;
    FN = 0;
    totP = 0;
    totN = 0;
    for j = 1:bns
        for k = j:bns
            jk = dags{i}(j,k);  kj = dags{i}(k,j);
            jk0 = dag0(j,k); kj0 = dag0(k,j);
            if jk + kj == 0
                totN = totN+1;
            else
                totP = totP+1;
            end
            switch jk0
                case 0
                    switch kj0
                        case 0
                            if ~(jk == 0 && kj == 0)
                                FP = FP + 1;
                            else
                                TN = TN + 1;
                            end
                        case 1
                            if     jk == 0 && kj == 0
                                FN = FN + 1;
                            elseif jk == 0 && kj == 1
                                TP = TP + 1;
                            elseif jk == 1 && kj == 0
                                TP2 = TP2 + 1;
                            elseif jk == 2 && kj == 2
                                TP2 = TP2 + 1;
                            end
                    end
                case 1
                    if     jk == 0 && kj == 0
                        FN = FN + 1;
                    elseif jk == 0 && kj == 1
                        TP2 = TP2 + 1;
                    elseif jk == 1 && kj == 0
                        TP = TP + 1;
                    elseif jk == 2 && kj == 2
                        TP2 = TP2 + 1;
                    end
            end
        end
    end
    if totP0==0
        se(i) = 1;
    else
        se(i) = (TP + TP2/2)/(TP + TP2/2 + FN);
    end
    if totP == 0
        pr(i) = 1;
    else
        pr(i) = (TP + TP2/2)/(TP + TP2/2 + FP);
    end
    if totN0 == 0
        sp(i) = 1;
    else
        sp(i) = TN/totN0;
    end
    if pr(i)+se(i) == 0
        f1(i) = 0;
    else
        f1(i) = 2*pr(i)*se(i)/(pr(i)+se(i));
    end
    shd(i) = FP + FN + TP2/2;
end
end
