function [ss] = get_CI_test(data,bnet,tol)
n = size(bnet.dag,1);
ss = xor(true(n),diag(true(1,n)));
for i = 1:n-1
    for j = i+1:n
        ci = cond_indep_chisquare(i,j,[],data,'LRT',tol,bnet.node_sizes);
        if ci == 1
            ss(i,j) = false; ss(j,i) = false;
        end
    end
end
end
