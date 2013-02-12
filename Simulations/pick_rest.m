function M2 = pick_rest(M,M1)

% This function finds the complement of a M1 where M is the full set

M2 = M;
for i=1: length(M1)
    M2(M2==M1(i)) = [];
end
