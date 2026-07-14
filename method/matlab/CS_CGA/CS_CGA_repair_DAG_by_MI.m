function [pop] = CS_CGA_repair_DAG_by_MI(pop,MI,MP)
N = size(pop,2);
n = size(pop{1},1);
for i = 1:N
    [loop,is_dag] = get_loop(pop{i});
    while ~is_dag
        [e0,e1] = find(loop);
        loop_size = size(e0,1);
        edge_MI = zeros(loop_size,1);
        for j = 1:loop_size
            edge_MI(j) = MI(e0(j),e1(j));
        end
        [~,index] = sort(edge_MI);
        ii = e0(index(1));     jj = e1(index(1));
        pop{i}(ii,jj) = false;
        [loop,is_dag] = get_loop(pop{i});
    end
    pop{i} = naive_limit_parents(MP,pop{i});
end
end
