function [g_best,g_best_score] = get_best(score,pop)
j_star = find(score==max(score),1);
g_best_score = score(j_star);
g_best = pop{j_star};
end
