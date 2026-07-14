function [l_val] = get_allele(e1,e2)
if ~e1 && ~e2
    l_val = 1;
elseif ~e1 && e2
    l_val = 2;
else
    l_val = 3;
end
end
