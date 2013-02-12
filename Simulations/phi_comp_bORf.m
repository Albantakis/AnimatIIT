function [phi_MIP prob prob_prod_MIP MIP network] = phi_comp_bORf(numerator,denom,whole_sys_state,network,bf)
% Larissa: for smart purviews, op_context is assumed 0, op_min is assumed
% bf is back/forward flag (back = 1, forward = 2)

op_normalize = network.options(6);
op_small_phi = network.options(4);

% global BRs, global FRs
% global BRs_check FRs_check
% global BRs_check2 FRs_check2
% eps = 1e-10;

num_nodes_denom = length(denom);
num_nodes_numerator = length(numerator);
%% unpartitioned transition repertoire
% prob_w_old = Rs{convi(numerator),convi(denom)};

% current = convi(numerator); other = convi(denom);
current = sum(2.^(numerator-1))+1; other = sum(2.^(denom-1))+1; %Larissa: This is probably a step where only binary works!

if (bf == 1)
    if isempty(network.BRs{current,other})
        network.BRs{current,other} = comp_pers_cpt(network.nodes,numerator,denom,whole_sys_state,'backward');
%             BRs_check{current,other} = comp_pers_single(numerator,denom,whole_sys_state,p,1);

%     if ~all(abs(BRs{current,other} - BRs_check{current,other}) <= eps)
%         disp('BR CHECK:')
%         disp(numerator)
%         disp(denom)
%         disp(BRs{current,other})
%         disp(BRs_check{current,other})
%         disp(BRs{current,other}(:) == BRs_check{current,other}(:))
%     end
    end
    prob_w = network.BRs{current,other};
elseif (bf == 2)
    if isempty(network.FRs{current,other})
        network.FRs{current,other} = comp_pers_cpt(network.nodes,numerator,denom,whole_sys_state,'forward');
%                 FRs_check{current,other} = comp_pers_single(numerator,denom,whole_sys_state,p,2);

%     if ~all(abs(FRs{current,other} - FRs_check{current,other}) <= eps)
%         disp('FR CHECKx:')
%         disp(numerator)
%         disp(denom)
%         disp(FRs{current,other})
%         disp(FRs_check{current,other})
%         disp(FRs_check2{current,other})
%         disp(FRs{current,other}(:) == FRs_check{current,other}(:))
%     end
    end
    prob_w = network.FRs{current,other};
end
    
%     disp('CHECK NEW COMPUTATION:');
%     disp('prob_w_old:')
%     disp(prob_w_old);
%     disp('prob_w:')
%     disp(prob_w);
%     disp(prob_w_old == prob_w);

    
prob = cell(2,1);
prob{bf} = prob_w(:);
if bf == 1 %backward is calculated, forward is maxent
    forward_maxent_dist = comp_pers_cpt(network.nodes,[],denom,[],'forward');
    prob{2} = forward_maxent_dist;
elseif bf == 2
    uniform_dist = ones(size(prob{bf}))/length(prob{bf});
    prob{1} = uniform_dist;
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

for i = 1:num_denom_partitions % past or future
    denom_part1 = denom_partitions1{i};
    denom_part2 = denom_partitions2{i};
    for j=1: num_numerator_partitions % present
        numerator_part1 = numerator_partitions1{j};
        numerator_part2 = numerator_partitions2{j};
        Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2);

        if Norm ~= 0
%             current_1 = convi(numerator_part1); current_2 = convi(numerator_part2);
%             other_1 = convi(denom_part1); other_2 = convi(denom_part2);
            current_1 = sum(2.^(numerator_part1-1))+1;
            current_2 = sum(2.^(numerator_part2-1))+1;
            other_1 = sum(2.^(denom_part1-1))+1;
            other_2 = sum(2.^(denom_part2-1))+1;
            
            if (bf == 1)
                if isempty(network.BRs{current_1,other_1})
                    network.BRs{current_1,other_1} = comp_pers_cpt(network.nodes,numerator_part1,denom_part1,whole_sys_state,'backward');
%                             BRs_check{current_1,other_1} = comp_pers_single(numerator_part1,denom_part1,whole_sys_state,p,1);

%     if ~all(abs(BRs{current_1,other_1} - BRs_check{current_1,other_1}) <= eps)
%         disp('BR CHECK:')
%         disp(numerator_part1)
%         disp(denom_part1)
%         disp(BRs{current_1,other_1})
%         disp(BRs_check{current_1,other_1})
%         disp(BRs{current_1,other_1}(:) == BRs_check{current_1,other_1}(:))
%     end
                end
                prob_p1 = network.BRs{current_1,other_1};
                
                if isempty(network.BRs{current_2,other_2})
                    network.BRs{current_2,other_2} = comp_pers_cpt(network.nodes,numerator_part2,denom_part2,whole_sys_state,'backward');
%                                                 BRs_check{current_2,other_2} = comp_pers_single(numerator_part2,denom_part2,whole_sys_state,p,1);

%     if ~all(abs(BRs{current_2,other_2} - BRs_check{current_2,other_2}) <= eps)
%         disp('BR CHECK:')
%         disp(numerator_part2)
%         disp(denom_part2)
%         disp(BRs{current_2,other_2})
%         disp(BRs_check{current_2,other_2})
%         disp(BRs{current_2,other_2}(:) == BRs_check{current_2,other_2}(:))
%     end
                end
                prob_p2 = network.BRs{current_2,other_2};
                
            elseif (bf == 2)
                if isempty(network.FRs{current_1,other_1})
                    network.FRs{current_1,other_1} = comp_pers_cpt(network.nodes,numerator_part1,denom_part1,whole_sys_state,'forward');
%                                             FRs_check{current_1,other_1} = comp_pers_single(numerator_part1,denom_part1,whole_sys_state,p,2);

%     if ~all(abs(FRs{current_1,other_1} - FRs_check{current_1,other_1}) <= eps)
%         disp('FR CHECK:')
%         disp(numerator_part1)
%         disp(denom_part1)
%         disp(FRs{current_1,other_1})
%         disp(FRs_check{current_1,other_1})
%         disp(FRs_check2{current_1,other_1})
%         disp(FRs{current_1,other_1}(:) == FRs_check{current_1,other_1}(:))
%     end
                end
                prob_p1 = network.FRs{current_1,other_1};
                
                if isempty(network.FRs{current_2,other_2})
                   network.FRs{current_2,other_2} = comp_pers_cpt(network.nodes,numerator_part2,denom_part2,whole_sys_state,'forward');
%                                                                     FRs_check{current_2,other_2} = comp_pers_single(numerator_part2,denom_part2,whole_sys_state,p,2);

%     if ~all(abs(FRs{current_2,other_2} - FRs_check{current_2,other_2}) <= eps)
%         disp('FR CHECK:')
%         disp(numerator_part2)
%         disp(denom_part2)
%         disp(FRs{current_2,other_2})
%         disp(FRs_check{current_2,other_2})
%         disp(FRs_check2{current_2,current_2})
%         disp(FRs{current_2,other_2}(:) == FRs_check{current_2,other_2}(:))
%     end
                end
                prob_p2 = network.FRs{current_2,other_2};
            end
            
%             prob_p = prob_prod_comp(prob_p1(:),prob_p2(:),denom,denom_part1,0);

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
            elseif (op_small_phi == 1)
                phi = L1norm(prob{bf},prob_p);
            elseif (op_small_phi == 2)
                phi = emd_hat_gd_metric_mex(prob{bf},prob_p,gen_dist_matrix(length(prob_p)));
%             elseif (op_small_phi == 3)
%                 phi = k_distance(prob{bf},prob_p);
            end
            
            
        else
            prob_prod_vec{i,j,bf} = [];
            phi = Inf;
        end
        
        if phi == 0
            phi_MIP = [0 0];
            prob_prod_MIP = cell(2,1);
            MIP = cell(2,2,2);
            return
        end

        phi_cand(i,j,bf,1) = phi;
        phi_cand(i,j,bf,2) = phi/Norm;

%             fprintf('phi=%f phi_norm=%f %s-%s -%s\n',phi,phi/Norm,mod_mat2str(xp_1),mod_mat2str(numerator_part1),mod_mat2str(xf_1));
    end
end

MIP = cell(2,2,2);
phi_MIP = zeros(1,2);
prob_prod_MIP = cell(2,1);

[phi_MIP(bf) i j] = min2(phi_cand(:,:,bf,1),phi_cand(:,:,bf,2),op_normalize);
prob_prod_MIP{bf} = prob_prod_vec{i,j,bf};

MIP{1,1,bf} = denom_partitions1{i};
MIP{2,1,bf} = denom_partitions2{i};
MIP{1,2,bf} = numerator_partitions1{j};
MIP{2,2,bf} = numerator_partitions2{j};

end

function Norm = Normalization(xp_1,xp_2,numerator_part1,numerator_part2,xf_1,xf_2)

if nargin == 4
    Norm = min(length(numerator_part1),length(xp_2)) + min(length(numerator_part2),length(xp_1));
else
    Norm = min(length(numerator_part1),length(xp_2)) + min(length(numerator_part2),length(xp_1)) ...
        + min(length(numerator_part1),length(xf_2)) + min(length(numerator_part2),length(xf_1));
end

end

function [X_min i_min j_min k_min] = min3(X,X2)
X_min = Inf; % minimum of normalized phi
X_min2 = Inf; % minimum of phi
i_min = 1;
j_min = 1;
k_min = 1;

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
            if dif > epsilon || abs(dif) < epsilon
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

