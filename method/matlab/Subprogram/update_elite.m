function [g_best,g_best_score,p,score] = update_elite(g_best,g_best_score,p,score)
[g_best2,g_best_score2] = get_best(score,p);
if g_best_score2 > g_best_score
    g_best = g_best2;
    g_best_score = g_best_score2;
else
    j_worst = find(score==min(score),1);
    p{j_worst} = g_best;
    score(j_worst) = g_best_score;
end
end
