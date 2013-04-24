function [phi_MIP, prob, prob_prod_MIP, MIP, network] = phi_comp_bf(subsystem,numerator,denom_past,denom_future,whole_sys_state,network)
% compute small phi of a given purview...?

op_normalize = network.options(6);
op_small_phi = network.options(4);
op_parfor = network.options(9);
op_extNodes = network.options(11);
op_complex = network.options(3);

num_nodes_denom_past = length(denom_past);
num_nodes_numerator = length(numerator);
num_nodes_denom_future = length(denom_future);

if op_parfor == 2 && op_extNodes == 0 && op_complex == 1
    BRs = network.BRs{subsystem2index(subsystem)};
    FRs = network.FRs{subsystem2index(subsystem)};
else    
    BRs = network.BRs;
    FRs = network.FRs;
end

if op_extNodes == 0
    extNodes = setdiff(network.full_system, subsystem);
else
    extNodes = [];
end   

%% unpartitioned transition repertoire
current = sum(2.^(numerator-1))+1; past = sum(2.^(denom_past-1))+1; future = sum(2.^(denom_future-1))+1;

if isempty(BRs{current,past})
    BRs{current,past} = comp_pers_cpt(network.nodes,numerator,denom_past,whole_sys_state,'backward',extNodes,network.past_state);
end
prob_bw = BRs{current,past};

if isempty(FRs{current,future})
    FRs{current,future} = comp_pers_cpt(network.nodes,numerator,denom_future,whole_sys_state,'forward',extNodes);
end
prob_fw = FRs{current,future};


prob = cell(2,1);
prob{1} = prob_bw(:);
prob{2} = prob_fw(:);

%% more than one
if num_nodes_denom_past ~= 0
    [denom_past_partitions_1 denom_past_partitions_2 num_denom_partitions] = bipartition(denom_past,num_nodes_denom_past); % partition of denom_past
else
    denom_past_partitions_1{1} = []; denom_past_partitions_2{1} = []; num_denom_partitions = 1;
end
[num_numerator_partitions1 num_numerator_partitions2 num_numerator_partitions] = bipartition(numerator,num_nodes_numerator,1); % partition of numerator

MIP = cell(2,2,2);
phi_MIP = zeros(1,2);
prob_prod_MIP = cell(2,1);

phi_cand = zeros(num_denom_partitions,num_numerator_partitions,2,2);
prob_prod_vec = cell(num_denom_partitions,num_numerator_partitions,2);

for bf = 1:2 % past and future
  phi_zero_found = 0;  
for i = 1:num_denom_partitions % past or future
    denom_part1 = denom_past_partitions_1{i};
    denom_part2 = denom_past_partitions_2{i};

    for j = 1:num_numerator_partitions % present
        numerator_part1 = num_numerator_partitions1{j};
        numerator_part2 = num_numerator_partitions2{j};
        
        Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2);
       
        current_1 = sum(2.^(numerator_part1-1))+1;
        current_2 = sum(2.^(numerator_part2-1))+1;
        other_1 = sum(2.^(denom_part1-1))+1;
        other_2 = sum(2.^(denom_part2-1))+1;

            if Norm ~= 0

                if bf == 1
                    if isempty(BRs{current_1,other_1})
                            BRs{current_1,other_1} = comp_pers_cpt(network.nodes,numerator_part1,denom_part1,whole_sys_state,'backward',extNodes,network.past_state);
                    end
                    prob_p1 = BRs{current_1,other_1};

                    if isempty(BRs{current_2,other_2})
                            BRs{current_2,other_2} = comp_pers_cpt(network.nodes,numerator_part2,denom_part2,whole_sys_state,'backward',extNodes,network.past_state);
                    end
                    prob_p2 = BRs{current_2,other_2};

                else

                    if isempty(FRs{current_1,other_1})
                        FRs{current_1,other_1} = comp_pers_cpt(network.nodes,numerator_part1,denom_part1,whole_sys_state,'forward',extNodes);
                    end
                    prob_p1 = FRs{current_1,other_1};

                    if isempty(FRs{current_2,other_2})
                        FRs{current_2,other_2} = comp_pers_cpt(network.nodes,numerator_part2,denom_part2,whole_sys_state,'forward',extNodes);
                    end
                    prob_p2 = FRs{current_2,other_2};

                end
                
                    if isempty(prob_p1)
                        prob_p = prob_p2(:);
                    elseif isempty(prob_p2)
                        prob_p = prob_p1(:);
                    else
                        prob_p_test = bsxfun(@times,prob_p1,prob_p2);
                        prob_p = prob_p_test(:);
                    end
                prob_prod_vec{i,j,bf} = prob_p;
                
                if (op_small_phi == 0)
                    phi = KLD(prob{bf},prob_p);
                elseif op_small_phi == 1
                    phi = L1norm(prob{bf},prob_p); 
                elseif op_small_phi == 2
                    %phi = emd_hat_gd_metric_mex(prob{bf},prob_p,gen_dist_matrix(length(prob_p)));                  
                    phi = emd_hat_gd_metric_mex(prob{bf},prob_p,network.gen_dist_matrix(1:length(prob_p),1:length(prob_p)));                  
                end
                
            else
                prob_prod_vec{i,j,bf} = [];
                phi = Inf;
            end
            
            if phi == 0
                phi_zero_found = 1;
                break
            end
            
            phi_cand(i,j,bf,1) = phi;
            phi_cand(i,j,bf,2) = phi/Norm;
    end
        if phi_zero_found
            break
        end
end
    
        if phi_zero_found
            phi_MIP(bf)  = 0;

        else 
            [phi_MIP(bf) i j] = min2(phi_cand(:,:,bf,1),phi_cand(:,:,bf,2),op_normalize);
            prob_prod_MIP{bf} = prob_prod_vec{i,j,bf};

            MIP{1,1,bf} = denom_past_partitions_1{i};
            MIP{2,1,bf} = denom_past_partitions_2{i};
            MIP{1,2,bf} = num_numerator_partitions1{j};
            MIP{2,2,bf} = num_numerator_partitions2{j};
        end
end

if op_parfor == 2 && op_extNodes == 0
    network.BRs{subsystem2index(subsystem)} = BRs;
    network.FRs{subsystem2index(subsystem)} = FRs;
else    
    network.BRs = BRs;
    network.FRs = FRs;
end
end

function [X_min i_min j_min k_min] = min3(X,X2,op_normalize)
X_min = Inf; % minimum of normalized phi (or unnormalized if op_normalize == 0)
X_min2 = Inf; % minimum of phi
i_min = 1;
j_min = 1;
k_min = 1;

if (op_normalize == 1)
    for i=1: size(X,1)
        for j=1: size(X,2)
            for k=1: size(X,3)
                if X(i,j,k) <= X_min && X2(i,j,k) <= X_min2
                    X_min = X(i,j,k);
                    X_min2 = X2(i,j,k);
                    i_min = i;
                    j_min = j;
                    k_min = k;
                end            
            end
        end
    end
else
    for i=1: size(X,1)
        for j=1: size(X,2)
            for k=1: size(X,3)
                if X2(i,j,k) <= X_min
    %                 X_min = X(i,j,k);
                    X_min = X2(i,j,k);
                    i_min = i;
                    j_min = j;
                    k_min = k;
                end            
            end
        end
    end
end
end


function [phi_min_choice, i_min, j_min] = min2(phi,phi_norm,op_normalize)
phi_norm_min = Inf; % minimum of normalized phi
phi_min = Inf; % minimum of phi
i_min = 1;
j_min = 1;
epsilon = 10^-10;

if (op_normalize == 1 || op_normalize == 2)
    for i=1: size(phi,1)
        for j=1: size(phi,2)
%             if phi_norm(i,j) <= phi_norm_min && phi(i,j) <= phi_min
            dif = phi_norm_min - phi_norm(i,j); 
            if dif > epsilon || abs(dif) < epsilon  %Larissa: instead of phi <= phi_min
                phi_min = phi(i,j);
                phi_norm_min = phi_norm(i,j);
                i_min = i;
                j_min = j;
            end
        end
    end
else
    for i=1: size(phi,1)
        for j=1: size(phi,2)
            dif = phi_min - phi(i,j); 
            if dif > epsilon || abs(dif) < epsilon
                phi_min = phi(i,j);
                phi_norm_min = phi_norm(i,j);
                i_min = i;
                j_min = j;
            end
        end
    end
end

if (op_normalize == 0 || op_normalize == 1)
    phi_min_choice = phi_min;
else
    phi_min_choice = phi_norm_min;
end

end
