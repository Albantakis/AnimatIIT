function [B1 B2 N_b] = bipartition(X,N_max,op)
% X: subset, B1: group #1, B2: group #2

if nargin < 3
    op = 0;
end

if op == 1
    N_max = floor(N_max/2);
end

N = length(X); % number of elements

N_b = 0; % number of bipartition
for i=0: N_max
    N_b = N_b + nchoosek(N,i);
end

B1 = cell(N_b,1);
B2 = cell(N_b,1);

i_b = 1;
for i=0: N_max
    if i== 0
        B1{i_b} = [];
        B2{i_b} = X;
        i_b = i_b + 1;
    else
        C_b = nchoosek(1:N,i);
        N_C = size(C_b,1); % this equals nchoosek(N,i)
        for j= 1:N_C
            B1{i_b} = X(C_b(j,:));
            B2_temp = 1:N;
            B2_temp(C_b(j,:)) = [];
            B2{i_b} = X(B2_temp);
            i_b = i_b + 1;
        end
    end
    
end