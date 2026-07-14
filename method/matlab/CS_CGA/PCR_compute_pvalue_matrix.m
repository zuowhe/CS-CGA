function [p_value] = PCR_compute_pvalue_matrix(data, bnet, tol)
    n = size(bnet.dag, 1);
    ss = xor(true(n), diag(true(1, n)));
    p_value = double(ss);
    for i = 1:n-1
        for j = i+1:n
            [~, ~, alpha2] = PCR_conditional_independence_pvalue(i, j, [], data, 'LRT', tol, bnet.node_sizes);
            p_value(i, j) = alpha2;  p_value(j, i) = alpha2;
        end
    end
end
