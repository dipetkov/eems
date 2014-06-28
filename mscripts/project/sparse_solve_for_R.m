

function [R,time] = sparse_solve_for_R(Mstruct,mValues)


tic( );
G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
R = resistance_distance(G);
time = toc( );

%% This is the McRae's approximation to ETab
%% when the coalescent scale is Na (as in ms)
%% rather than N = sum_a Na

R = R/4 + 1;
