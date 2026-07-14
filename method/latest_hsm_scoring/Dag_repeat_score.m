function [score,pop_history, score_history,repeated_count1,length_history,cache] = Dag_repeat_score(data, N2, ns, pop, pop_history, score_history, length_history,scoring_fn, cache)
% 检查当前代的个体是否已经计算过评分，避免重复计算
    score = zeros(1, N2);  % 初始化当前代的评分
    repeated_count1 = 0;
    for j = 1:N2
        repeated = false;
        % 对每个个体检查是否已经计算过评分
        for k = 1:length_history
            if isequal(pop{j}, pop_history{k})  % 如果当前个体在历史种群中已经存在
                score(j) = score_history(k);  % 使用已保存的评分
                repeated = true;
                repeated_count1 = repeated_count1 + 1;
                break;
            end
        end
        
        if ~repeated
            % 计算该个体的评分并保存
            [individual_score, cache] = score_dags(data, ns, pop(j), 'scoring_fn', scoring_fn, 'cache', cache);
            score(j) = individual_score;
            % 保存当前个体和其评分
            length_history = length_history + 1;
            pop_history{length_history} = pop{j};
            score_history(length_history) = individual_score;
        end
    end
end

