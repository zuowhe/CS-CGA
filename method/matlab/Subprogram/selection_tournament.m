function [pop1,score1] = selection_tournament(N,N2,pop,score,tour)
i = 0;
pop1 = cell(1,N2);
score1 = -Inf * ones(1,N2);
while i < N
    index = ceil(rand(1,tour) *N2);
    score_list = score(index);
    best_list = index(find(score_list == max(score_list)));
    best_pop = best_list(1);
    i = i+1;
    pop1{i} = pop{best_pop};
    score1(i) = score(best_pop);
end
end
