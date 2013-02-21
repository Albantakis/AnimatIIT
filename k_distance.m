function [k_diff, intra_1, intra_2, inter] = k_distance(P,Q)

%% compute Intra and interdistance based on the Kernel idea, for distributions
% P,Q : Two distributions


%% 

% if (sum(P)~= 1) || (sum(Q)~= 1)
%     warning('Either or both distributions do not add up to 1');
% end 

if min(min(P),min(Q)) < 0
    warning('Either or both distributions contain negative values');
end 

if (length(P)~= length(Q))
    warning('Distributions do not have the same support');
end

%compute interdistance:
inter = k_inter(P,P) + k_inter(Q,Q) - 2* k_inter(P,Q);

%compute intradistances:
intra_1 = k_intra(P);
intra_2 = k_intra(Q);

k_diff = (intra_1 + intra_2) * inter;

end