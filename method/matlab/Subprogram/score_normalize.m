function [norm_score,sum_norm_score] = score_normalize(score,sc_star,compute_sum)
eps = 0.0000001;
sum_norm_score = 0;
sc_worst = min(score);
range = sc_star-sc_worst;
norm_score = (score-sc_worst)/(eps+range);
if compute_sum
    sum_norm_score = sum(norm_score)+eps;
end
end
