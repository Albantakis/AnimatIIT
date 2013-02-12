function [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_M, phi_M,MIP_M,subset,M_i_max)
%% Irreducible points
if M_i_max
    REP_cell = prob_M{M_i_max,1};
    REP_prod_cell = prob_M{M_i_max,2};
    MIP_cell = MIP_M{M_i_max};
    IRR_phi = phi_M{M_i_max}(:,1);
else
    REP_cell = prob_M{1};
    REP_prod_cell = prob_M{2};
    MIP_cell = MIP_M;
    IRR_phi = phi_M(:,1);   
end

N = length(subset);

index_vec_IRR = find(IRR_phi ~= 0);
N_IRR = length(index_vec_IRR);


IRR_REP = cell(N_IRR,2);
for i=1: N_IRR
    j = index_vec_IRR(i);
    IRR_REP{i,1} = REP_cell{j};
    IRR_REP{i,2} = REP_prod_cell{j};
end


IRR_MIP = cell(N_IRR,1);
for i=1: N_IRR
    j = index_vec_IRR(i);
    IRR_MIP{i} = MIP_cell{j};
end

% IRR_REP = prob_M{M_i_max};
% IRR_REP(:,IRR_phi==0) = [];

% IRR_phi(IRR_phi(:,1) == 0) = [];
IRR_phi = phi_M{M_i_max}(index_vec_IRR,:);


M_cell = cell(2^N-1,1);
k = 1;
for i=1: N
    C = nchoosek(subset,i);
    N_C = size(C,1);
    for j=1: N_C
        C_j = C(j,:);   
        M_cell{k} = C_j;
        k = k + 1;
    end
end

M_IRR = cell(N_IRR,1);

for i=1: N_IRR
    j = index_vec_IRR(i);
    M_IRR{i} = M_cell{j};
end