function [loop,is_dag] = get_loop(p)
n = size(p,1);
loop = p;
change = true;
while change
    change = false;
    for i = 1:n
        if ~any(loop(i,:)) && any(loop(:,i))
            loop(:,i) = 0;
            change = true;
        end
        if ~any(loop(:,i)) && any(loop(i,:))
            loop(i,:) = 0;
            change = true;
        end
    end
end
if any(loop,'all')
    is_dag = false;
else
    is_dag = true;
end
end
