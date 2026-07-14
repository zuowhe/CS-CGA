function score_gap = score_gap_normalize(current_best_score, historical_best_score, population_scores)
scores = population_scores(:);
score_range = max(scores) - min(scores);
eps_scale = 1e-12 * max(1, max(abs(scores)));
score_scale = max(score_range, eps_scale);
score_gap = (current_best_score - historical_best_score) / score_scale;
end
