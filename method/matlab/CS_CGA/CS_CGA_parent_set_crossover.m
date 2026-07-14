function [s] = CS_CGA_parent_set_crossover(N,s)
n = size(s{1},1);
for i=1:N
    j = i;
    while i==j
        j = randi(N);
    end
    s{N+i} = s{i};
    for ps=1:n
        if round(rand)
            s{N+i}(:,ps) = s{j}(:,ps);
        end
    end
end
end
