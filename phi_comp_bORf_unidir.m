function [phi_MIP prob network] = phi_comp_bORf_unidir(M1,M2,numerator,denom,whole_sys_state,network,bf_option,bfcut_option)
% compute small phi of a given higher order purview...?
% if the system is cut unidirectionally

% options = the options
% subsystem = a system
% numerator = state of the system??
% denom_past = 
% denom_future = 
% whole_sys_state = 
% p = TPM as a 2^N x N matrix
% b_table
% BRs
% FRs
bf = 1; %Larissa: That's only there to keep the dimensions as before, should be updated!

op_normalize = network.options(6);
op_small_phi = network.options(4);


num_nodes_denom = length(denom);
num_nodes_numerator = length(numerator);

%% unpartitioned transition repertoire

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rule: M1 <- M2 noised --> past: M1M2/M1 x M2/M1M2 future: M1/M1M2 x M1M2/M2%
%      M1 -> M2 noised --> future: M1/M1M2 X M1M1/M2 past: M1M2/M1 x M2/M1M2%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%M1 <- M2 noised and past, or M1 -> M2 noised and future
if (strcmp(bfcut_option,'BRcut') && strcmp(bf_option,'backward')) || (strcmp(bfcut_option,'FRcut') && strcmp(bf_option,'forward')) 
    numeratorM1 = numerator;
    numeratorM2 = numerator(ismember(numerator, M2));
    denomM1 = denom(ismember(denom,M1));
    denomM2 = denom;
%M1 <- M2 noised and future, or M1 -> M2 noised and past
elseif (strcmp(bfcut_option,'BRcut') && strcmp(bf_option,'forward')) || (strcmp(bfcut_option,'FRcut') && strcmp(bf_option,'backward'))   
    numeratorM1 = numerator(ismember(numerator, M1));
    numeratorM2 = numerator;
    denomM1 = denom;
    denomM2 = denom(ismember(denom,M2));
end
currentM1 = sum(2.^(numeratorM1-1))+1;  %index of BR/FR matrix
currentM2 = sum(2.^(numeratorM2-1))+1;  %will be 1 if empty!
otherM1 = sum(2.^(denomM1-1))+1;    
otherM2 = sum(2.^(denomM2-1))+1;    

if strcmp(bf_option,'backward')
    if isempty(network.BRs{currentM1,otherM1}) && otherM1 > 1
        network.BRs{currentM1,otherM1} = comp_pers_cpt(network.nodes,numeratorM1,denomM1,whole_sys_state,bf_option);
    end
    prob_M1 = network.BRs{currentM1,otherM1};
    
    if isempty(network.BRs{currentM2,otherM2}) && otherM2 > 1
       network.BRs{currentM2,otherM2} = comp_pers_cpt(network.nodes,numeratorM2,denomM2,whole_sys_state,bf_option);
    end
    prob_M2 = network.BRs{currentM2,otherM2};
elseif strcmp(bf_option,'forward')
    if isempty(network.FRs{currentM1,otherM1}) && otherM1 > 1
        network.FRs{currentM1,otherM1} = comp_pers_cpt(network.nodes,numeratorM1,denomM1,whole_sys_state,bf_option);
    end
    prob_M1 = network.FRs{currentM1,otherM1};
   
    if isempty(network.FRs{currentM2,otherM2}) && otherM2 > 1
       network.FRs{currentM2,otherM2} = comp_pers_cpt(network.nodes,numeratorM2,denomM2,whole_sys_state,bf_option);
    end
    prob_M2 = network.FRs{currentM2,otherM2};
end    
   
if isempty(prob_M1)
    prob = prob_M2(:);
elseif isempty(prob_M2)
    prob = prob_M1(:);
else
    prob_test = bsxfun(@times,prob_M1,prob_M2);
    prob_test = prob_test ./ sum(prob_test(:));
    prob = prob_test(:);
end 
%% more than one
if num_nodes_denom ~= 0
    [denom_partitions1 denom_partitions2 num_denom_partitions] = bipartition(denom,num_nodes_denom); % partition of xp
else
    denom_partitions1{1} = []; denom_partitions2{1} = []; num_denom_partitions = 1;
end

[numerator_partitions1 numerator_partitions2 num_numerator_partitions] = bipartition(numerator,num_nodes_numerator,1); % partition of numerator

phi_cand = zeros(num_denom_partitions,num_numerator_partitions,2,2);
prob_prod_vec = cell(num_denom_partitions,num_numerator_partitions,2,2);

phi_zero_found = 0;  
for i = 1:num_denom_partitions % past or future
    denom_part1 = denom_partitions1{i};
    denom_part2 = denom_partitions2{i};
    
    for j=1: num_numerator_partitions % present
        numerator_part1 = numerator_partitions1{j};
        numerator_part2 = numerator_partitions2{j};
               
        Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2);
        
        if Norm ~= 0
            
            %M1 <- M2 noised and past, or M1 -> M2 noised and future
            if (strcmp(bfcut_option,'BRcut') && strcmp(bf_option,'backward')) || (strcmp(bfcut_option,'FRcut') && strcmp(bf_option,'forward')) 
                numeratorM1_part1 = numerator_part1;
                numeratorM2_part1 = numerator_part1(ismember(numerator_part1, M2));
                numeratorM1_part2 = numerator_part2;
                numeratorM2_part2 = numerator_part2(ismember(numerator_part2, M2));
                denomM1_part1 = denom_part1(ismember(denom_part1,M1));
                denomM2_part1 = denom_part1;
                denomM1_part2 = denom_part2(ismember(denom_part2,M1));
                denomM2_part2 = denom_part2;
            %M1 <- M2 noised and future, or M1 -> M2 noised and past
            elseif (strcmp(bfcut_option,'BRcut') && strcmp(bf_option,'forward')) || (strcmp(bfcut_option,'FRcut') && strcmp(bf_option,'backward'))   
                numeratorM1_part1 = numerator_part1(ismember(numerator_part1, M1));
                numeratorM2_part1 = numerator_part1;
                numeratorM1_part2 = numerator_part2(ismember(numerator_part2, M1));
                numeratorM2_part2 = numerator_part2;
                denomM1_part1 = denom_part1;
                denomM2_part1 = denom_part1(ismember(denom_part1,M2));
                denomM1_part2 = denom_part2;
                denomM2_part2 = denom_part2(ismember(denom_part2,M2));
            end
            
            currentM1_1 = sum(2.^(numeratorM1_part1-1))+1;
            currentM2_1 = sum(2.^(numeratorM2_part1-1))+1;
            currentM1_2 = sum(2.^(numeratorM1_part2-1))+1;
            currentM2_2 = sum(2.^(numeratorM2_part2-1))+1;
            
            otherM1_1 = sum(2.^(denomM1_part1-1))+1;
            otherM2_1 = sum(2.^(denomM2_part1-1))+1;
            otherM1_2 = sum(2.^(denomM1_part2-1))+1;
            otherM2_2 = sum(2.^(denomM2_part2-1))+1;
        
            if strcmp(bf_option,'backward')
                if isempty(network.BRs{currentM1_1,otherM1_1}) && otherM1_1 > 1
                    network.BRs{currentM1_1,otherM1_1} = comp_pers_cpt(network.nodes,numeratorM1_part1,denomM1_part1,whole_sys_state,bf_option);
                end
                probM1_p1 = network.BRs{currentM1_1,otherM1_1};
                
                if isempty(network.BRs{currentM1_2,otherM1_2}) && otherM1_2 > 1
                    network.BRs{currentM1_2,otherM1_2} = comp_pers_cpt(network.nodes,numeratorM1_part2,denomM1_part2,whole_sys_state,bf_option);
                end
                probM1_p2 = network.BRs{currentM1_2,otherM1_2};
                
                if isempty(network.BRs{currentM2_1,otherM2_1}) && otherM2_1 > 1
                    network.BRs{currentM2_1,otherM2_1} = comp_pers_cpt(network.nodes,numeratorM2_part1,denomM2_part1,whole_sys_state,bf_option);
                end
                probM2_p1 = network.BRs{currentM2_1,otherM2_1};
                
                if isempty(network.BRs{currentM2_2,otherM2_2}) && otherM2_2 > 1
                    network.BRs{currentM2_2,otherM2_2} = comp_pers_cpt(network.nodes,numeratorM2_part2,denomM2_part2,whole_sys_state,bf_option);
                end
                probM2_p2 = network.BRs{currentM2_2,otherM2_2};

            elseif strcmp(bf_option,'forward')

                if isempty(network.FRs{currentM1_1,otherM1_1}) && otherM1_1 > 1
                    network.FRs{currentM1_1,otherM1_1} = comp_pers_cpt(network.nodes,numeratorM1_part1,denomM1_part1,whole_sys_state,bf_option);
                end
                probM1_p1 = network.FRs{currentM1_1,otherM1_1};
                
                if isempty(network.FRs{currentM1_2,otherM1_2}) && otherM1_2 > 1
                    network.FRs{currentM1_2,otherM1_2} = comp_pers_cpt(network.nodes,numeratorM1_part2,denomM1_part2,whole_sys_state,bf_option);
                end
                probM1_p2 = network.FRs{currentM1_2,otherM1_2};
                
                if isempty(network.FRs{currentM2_1,otherM2_1}) && otherM2_1 > 1
                    network.FRs{currentM2_1,otherM2_1} = comp_pers_cpt(network.nodes,numeratorM2_part1,denomM2_part1,whole_sys_state,bf_option);
                end
                probM2_p1 = network.FRs{currentM2_1,otherM2_1};
                
                if isempty(network.FRs{currentM2_2,otherM2_2}) && otherM2_2 > 1
                    network.FRs{currentM2_2,otherM2_2} = comp_pers_cpt(network.nodes,numeratorM2_part2,denomM2_part2,whole_sys_state,bf_option);
                end
                probM2_p2 = network.FRs{currentM2_2,otherM2_2};

            end
            
            % first find part 1
            %Don't flatten matrices here, because bsxfun below needs it in high dimension
            if isempty(probM1_p1)
                prob_p1 = probM2_p1;
            elseif isempty(probM2_p1)
                prob_p1 = probM1_p1;
            else
                prob_p1 = bsxfun(@times,probM1_p1,probM2_p1); 
                prob_p1 = prob_p1 ./ sum(prob_p1(:));
            end
            % then part 2
             %Don't flatten matrices here, because bsxfun below needs it in high dimension
            if isempty(probM1_p2)
                prob_p2 = probM2_p2;
            elseif isempty(probM2_p2)
                prob_p2 = probM1_p2;
            else
                prob_p2 = bsxfun(@times,probM1_p2,probM2_p2);
                prob_p2 = prob_p2 ./ sum(prob_p2(:));
            end
            % then multiply partition
            if isempty(prob_p1)
                prob_p = prob_p2(:);
            elseif isempty(prob_p2)
                prob_p = prob_p1(:);
            else
                prob_p_test = bsxfun(@times,prob_p1,prob_p2);
                prob_p_test = prob_p_test ./ sum(prob_p_test(:));
                prob_p = prob_p_test(:);
            end
            
            prob_prod_vec{i,j,bf} = prob_p;

            if (op_small_phi == 0)
                phi = KLD(prob,prob_p);
            elseif (op_small_phi == 1)
                phi = L1norm(prob,prob_p);
            elseif op_small_phi == 2
                phi = emd_hat_gd_metric_mex(prob,prob_p,gen_dist_matrix(length(prob_p)));
%             elseif (op_small_phi == 3)
%                 phi = k_distance(prob,prob_p);
            end

        else
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
    phi_MIP = 0;
else 
    [phi_MIP i j] = min2(phi_cand(:,:,bf,1),phi_cand(:,:,bf,2),op_normalize);
%     prob_prod_MIP{bf} = prob_prod_vec{i,j,bf};
% 
%     MIP{1,1,bf} = denom_past_partitions_1{i};
%     MIP{2,1,bf} = denom_past_partitions_2{i};
%     MIP{1,2,bf} = num_numerator_partitions1{j};
%     MIP{2,2,bf} = num_numerator_partitions2{j};
end
end

%% subfunctions

function Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2,xf_1,xf_2)

if nargin == 4
    Norm = min(length(numerator_part1),length(denom_part2)) + min(length(numerator_part2),length(denom_part1));
else
    Norm = min(length(numerator_part1),length(denom_part2)) + min(length(numerator_part2),length(denom_part1)) ...
        + min(length(numerator_part1),length(xf_2)) + min(length(numerator_part2),length(xf_1));
end

end

function [phi_min_choice i_min j_min] = min2(phi,phi_norm,op_normalize)
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
