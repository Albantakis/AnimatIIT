function M_i = trans_M(M,N)

% M is a subset

M_b = zeros(N,1);

M_l = length(M);

for i=1: M_l
    M_b(M(i)) = 1;
end

M_i = trans10(M_b)-1;