function [sc,conv] = CS_CGA_finalize_output(dag,data,bnet,scoring_fn,M,conv,last_gen)
[f1,se,sp] = eval_dags({dag},bnet.dag,1);
sc = score_dags(data,bnet.node_sizes,{dag},'scoring_fn',scoring_fn);
conv.f1(last_gen:M) = f1;
conv.se(last_gen:M) = se;
conv.sp(last_gen:M) = sp;
conv.sc(last_gen:M) = sc;
end
