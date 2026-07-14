function [CI_new, count_CIchange] = PCR_update_alpha(p_avg, ~, M, nodesNum, Dif_BIC, CI_new, count_CIchange)
CI_initial = p_avg / nodesNum;
if CI_new == 0
    CI_new = CI_initial;
    return;
end
if Dif_BIC < 0
    % PCR uses the normalized score gap to decide how strongly to relax alpha.
    eta = 1 / (1 + exp(Dif_BIC));
    if isnan(eta) || isinf(eta)
        eta = 0;
    end
    max_step = (1 - CI_initial) / max(M, 1);
    CI_new = CI_new + eta * max_step;
    CI_new = min(CI_new, 1);
    count_CIchange = count_CIchange + 1;
end
end
