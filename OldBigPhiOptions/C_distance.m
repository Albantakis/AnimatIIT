function [C_diff, intra_1, intra_2, inter] = C_distance(C,D)

%% compute Intra and interdistance based on the Kernel idea, for constellations
% C,D : Two constellations
% each constellation is a N X 2 cell array where each row is a concept and
% the first column is the distribution (point) and the second column is the
% phi value

% options:  

No_Nu = 0;   % a small amount to be added to avoid multiplying by zero

%% 


%compute interdistance:
inter = C_inter(C,C) + C_inter(D,D) - 2* C_inter(C,D);

%compute intradistances:
intra_1 = C_intra(C);
intra_2 = C_intra(D);

C_diff = (intra_1 + intra_2 + No_Nu) * inter;

end