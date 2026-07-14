function key = CS_CGA_dag_key(dag_matrix)
[row, col] = find(dag_matrix == 1);
key = generate_unique_key([row, col]);
end
function key = generate_unique_key(edges)
if isempty(edges)
    key = 'empty';
    return;
end
sorted_edges = sortrows(edges);
key_parts = cell(size(sorted_edges, 1), 1);
for i = 1:size(sorted_edges, 1)
    key_parts{i} = sprintf('%d_%d', sorted_edges(i,1), sorted_edges(i,2));
end
key = strjoin(key_parts, '_');
end
