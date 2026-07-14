function [x] = naive_limit_parents(MP,x)
EP = sum(x) - MP;
for v = 1:size(x,1)
    if EP(v) > 0
        x(:,v) = random_parent_limitation(EP(v),x(:,v));
    end
end
end
