function [max_time] = get_max_realtime(BN_NodesNum)
if BN_NodesNum<10
    max_time = 150;
elseif BN_NodesNum<30
    max_time = 5000;
elseif BN_NodesNum<40
    max_time = 5000;
elseif BN_NodesNum<50
    max_time = 8000;
elseif BN_NodesNum<60
    max_time = 30000;
elseif BN_NodesNum<80
    max_time = 30000;
elseif BN_NodesNum<120
    max_time = 40000;
else
    max_time = 40000;
end
end
