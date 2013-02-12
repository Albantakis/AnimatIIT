function [phi_MIP prob prob_prod_MIP MIP network] = phi_comp_bf(subsystem,numerator,denom_past,denom_future,whole_sys_state,network)
% compute small phi of a given purview...?
% 
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


op_normalize = network.options(6);
op_small_phi = network.options(4);

% global FRs, global b_table
% global BRs_check, global FRs_check
% global BRs_check2, global FRs_check2
% eps = 1e-10;
% global func_time, global inline_time
% global cpt_time tpm_time
% global nodes

num_nodes_denom_past = length(denom_past);
num_nodes_numerator = length(numerator);
num_nodes_denom_future = length(denom_future);

%% unpartitioned transition repertoire

% BRs = cell(2^N,2^N);
% current_convi = convi(numerator); past_convi = convi(denom_past); future_convi = convi(denom_future);
current = sum(2.^(numerator-1))+1; past = sum(2.^(denom_past-1))+1; future = sum(2.^(denom_future-1))+1;

% if (current_convi ~= current) || (past_convi ~= past) || (future_convi ~= future)
%     disp('ERROR!')
% end

% current = numerator; past = denom_past; future = denom_future;

if isempty(network.BRs{current,past})
%     tic
    network.BRs{current,past} = comp_pers_cpt(network.nodes,numerator,denom_past,whole_sys_state,'backward');
% cpt_time = cpt_time + toc;
% disp('old')
% tic
%     network.BRs_check{current,past} = comp_pers_single(numerator,denom_past,whole_sys_state,p,1);
% tpm_time = tpm_time + toc;
%     BRs{current,past} = comp_pers_single(numerator,denom_past,whole_sys_state,p,1);
    
%     if ~all(abs(BRs{current,past} - BRs_check{current,past}) <= eps)
%         disp('BR CHECK:')
%         disp(numerator)
%         disp(denom_past)
%         disp(BRs{current,past})
%         disp(BRs_check{current,past})
%         disp(BRs{current,past}(:) == BRs_check{current,past}(:))
%     end
end
prob_bw = network.BRs{current,past};

if isempty(network.FRs{current,future})

% tic
    network.FRs{current,future} = comp_pers_cpt(network.nodes,numerator,denom_future,whole_sys_state,'forward');
% cpt_time = cpt_time + toc;
% tic
%     FRs_check{current,future} = comp_pers_single(numerator,denom_future,whole_sys_state,p,2);
% tpm_time = tpm_time + toc;
    
%     disp('new result:')
%     disp(size(FRs{current,future}))
% %     disp(FRs{current,future}(:))
%     disp('check result:')
%     disp(size(FRs_check{current,future}))
%     disp(FRs_check{current,future})

%     if ~all(abs(FRs{current,future} - FRs_check{current,future}) <= eps)
%         disp('FR CHECK:')
%         disp(numerator)
%         disp(denom_past)
%         disp(FRs{current,future})
%         disp(FRs_check{current,future})
%         disp(FRs_check2{current,future})
%         disp(FRs{current,future}(:) == FRs_check{current,future}(:))
%     end
end
prob_fw = network.FRs{current,future};


prob = cell(2,1);
prob{1} = prob_bw(:);
prob{2} = prob_fw(:);

%% more than one
if num_nodes_denom_past ~= 0
    [denom_past_partitions_1 denom_past_partitions_2 num_denom_partitions] = bipartition(denom_past,num_nodes_denom_past); % partition of denom_past
else
    denom_past_partitions_1{1} = []; denom_past_partitions_2{1} = []; num_denom_partitions = 1;
end
% if num_nodes_denom_future ~= 0
%     [xf_b1 xf_b2 Nf_b] = bipartition(denom_future,num_nodes_denom_future); % partition of denom_future
% else
%     xf_b1{1} = []; xf_b2{1} = []; Nf_b = 1;
% end
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
        
%         if (all(numerator == [1 2]) && all(denom_past == [1 2]))
%             disp('numerator_part1:')
%             disp(numerator_part1)
%             disp('numerator_part2:')
%             disp(numerator_part2)
%             disp('denom_part1:')
%             disp(denom_part1)
%             disp('denom_part2:')
%             disp(denom_part2)
%             disp('Norm:')
%             disp(Norm)
%         end
%         current_1_convi = convi(numerator_part1);
%         current_2_convi = convi(numerator_part2);
%         other_1_convi = convi(denom_part1);
%         other_2_convi = convi(denom_part2);
        
        current_1 = sum(2.^(numerator_part1-1))+1;
        current_2 = sum(2.^(numerator_part2-1))+1;
        other_1 = sum(2.^(denom_part1-1))+1;
        other_2 = sum(2.^(denom_part2-1))+1;
        
%         if (current_1_convi ~= current_1) || (current_2_convi ~= current_2) || (other_1_convi ~= other_1) || (other_2 ~= other_2_convi)
%             disp('ERROR!')
%         end

        
%         current_1 = numerator_part1;
%         current_2 = numerator_part2;
%         other_1 = denom_part1;
%         other_2 = denom_part2;

        

            if Norm ~= 0

                if bf == 1
                    if isempty(network.BRs{current_1,other_1})
%                             tic
                            network.BRs{current_1,other_1} = comp_pers_cpt(network.nodes,numerator_part1,denom_part1,whole_sys_state,'backward');
%                             cpt_time = cpt_time + toc;
%                             % disp('old')
%                             tic
%                                 BRs_check{current_1,other_1} = comp_pers_single(numerator_part1,denom_part1,whole_sys_state,p,1);
%                             tpm_time = tpm_time + toc;
 
%                         if ~all(abs(BRs{current_1,other_1} - BRs_check{current_1,other_1}) <= eps)
%                         disp('BR CHECK:')
%                         disp(numerator_part1)
%                         disp(denom_part1)
%                         disp(BRs{current_1,other_1})
%                         disp(BRs_check{current_1,other_1})
%                         disp(BRs{current_1,other_1}(:) == BRs_check{current_1,other_1}(:))
%                         end
                    end
                    prob_p1 = network.BRs{current_1,other_1};

                    if isempty(network.BRs{current_2,other_2})
%                             tic
                            network.BRs{current_2,other_2} = comp_pers_cpt(network.nodes,numerator_part2,denom_part2,whole_sys_state,'backward');
%                             cpt_time = cpt_time + toc;
%                             % disp('old')
%                             tic
%                                 BRs_check{current_2,other_2} = comp_pers_single(numerator_part2,denom_part2,whole_sys_state,p,1);
%                             tpm_time = tpm_time + toc;
%                         if ~all(abs(BRs{current_2,other_2} - BRs_check{current_2,other_2}) <= eps)
%                             disp('BR CHECK:')
%                             disp(numerator_part2)
%                             disp(denom_part2)
%                             disp(BRs{current_2,other_2})
%                             disp(BRs_check{current_2,other_2})
%                             disp(BRs{current_2,other_2}(:) == BRs_check{current_2,other_2}(:))
%                         end
                    end
                    prob_p2 = network.BRs{current_2,other_2};

                else

                    if isempty(network.FRs{current_1,other_1})
%                         tic
                        network.FRs{current_1,other_1} = comp_pers_cpt(network.nodes,numerator_part1,denom_part1,whole_sys_state,'forward');
%                         cpt_time = cpt_time + toc;
%                         tic
%                         FRs_check{current_1,other_1} = comp_pers_single(numerator_part1,denom_part1,whole_sys_state,p,2);
%                         tpm_time = tpm_time + toc;
%     
%                         if ~all(abs(FRs{current_1,other_1} - FRs_check{current_1,other_1}) <= eps)
%                             disp('FR CHECK:')
%                             disp(numerator_part1)
%                             disp(denom_part1)
% %                             disp('new result')
%                             disp(FRs{current_1,other_1})
% %                             disp('old result')
%                             disp(FRs_check{current_1,other_1})
%                             disp(FRs_check2{current_1,other_1})
%                             disp(FRs{current_1,other_1}(:) == FRs_check{current_1,other_1}(:))
%                         end
                    end
                    prob_p1 = network.FRs{current_1,other_1};

                    if isempty(network.FRs{current_2,other_2})
%                         tic
                        network.FRs{current_2,other_2} = comp_pers_cpt(network.nodes,numerator_part2,denom_part2,whole_sys_state,'forward');
%                         cpt_time = cpt_time + toc;
%                         tic
%                         FRs_check{current_2,other_2} = comp_pers_single(numerator_part2,denom_part2,whole_sys_state,p,2);
%                         tpm_time = tpm_time + toc;
%     
%                         if ~all(abs(FRs{current_2,other_2} - FRs_check{current_2,other_2}) <= eps)    
%                             disp('FR CHECK:')
%                             disp(numerator_part2)
%                             disp(denom_part2)
%                             disp(FRs{current_2,other_2})
%                             disp(FRs_check{current_2,other_2})
%                             disp(FRs_check2{current_2,other_2})
%                             disp(FRs{current_2,other_2}(:) == FRs_check{current_2,other_2}(:))
%                         end
                    end
                    prob_p2 = network.FRs{current_2,other_2};

                end
                
%                 if exist('prob_prod_comp','file') == 3
% %                     tic
%                     prob_p = prob_prod_comp(prob_p1(:),prob_p2(:),denom_past,denom_part1,0); % ADDED (:)
% %                     func_time = func_time + toc;
%                 
%                 else    
%                     tic
                    if isempty(prob_p1)
                        prob_p = prob_p2(:);
                    elseif isempty(prob_p2)
                        prob_p = prob_p1(:);
                    else
    %                     prob_p_test = (expand_prob(prob_p1,denom_past,denom_part1) .* expand_prob(prob_p2,denom_past,denom_part2));
    %                     prob_p_test = prob_p_test(:)/sum(prob_p_test);
    %                     prob_p1_reshape = reshape(prob_p1,

    %                     repmat_vec = ones(1,N);
    %                     repmat_vec(denom_part2) = 2;
    %                     prob_p1_rep = repmat(prob_p1,repmat_vec);
    %                     repmat_vec = ones(1,N);
    %                     repmat_vec(denom_part1) = 2;
    %                     prob_p2_rep = repmat(prob_p2,repmat_vec) ;                 
    %                     prob_p_test = prob_p1_rep .* prob_p2_rep;
    %                     prob_p_test = prob_p_test(:);
                        prob_p_test = bsxfun(@times,prob_p1,prob_p2);
                        prob_p = prob_p_test(:);

    %                     if ~all(prob_p_test == prob_p_test2(:))
    %                         disp('***NOPE***')
    %                     end

    %                     prob_p_test = prob_p1(:) * prob_p2(:)';
    %                     prob_p_test = prob_p_test(:)
                    end
%                     cpt_time = cpt_time + toc;
%                 end
                
                
%                 if ~all(abs(prob_p - prob_p_test) <= 1e-5)
%                     disp('PROB CHECK:')
%                 disp(abs(prob_p - prob_p_test) <= 1e-5)
%                     disp('****************')
%                     disp('****************')
%                     disp('present part 1:')
%                     disp(numerator_part1)
%                     disp('other part 1:')
%                     disp(denom_part1)
%                     disp('present part 2:')
%                     disp(numerator_part2)
%                     disp('other part 2:')
%                     disp(denom_part2)
%                     disp(prob_p1)
%                     disp(prob_p2)
%                     disp(prob_p)
%                     disp(prob_p_test)
%                     disp('****************')
% %                 else
% %                     disp('present part 1:')
% %                     disp(numerator_part1)
% %                     disp('other part 1:')
% %                     disp(denom_part1)
% %                     disp('present part 2:')
% %                     disp(numerator_part2)
% %                     disp('other part 2:')
% %                     disp(denom_part2)
% %                     disp(prob_p1)
% %                     disp(prob_p2)
% %                     disp(prob_p)
% %                     disp(prob_p_test)
%                 end
                prob_prod_vec{i,j,bf} = prob_p;
                
                if (op_small_phi == 0)
                    phi = KLD(prob{bf},prob_p);
%                                     phi2 = KLD_old(prob{bf},prob_p);
%                 if (phi ~= phi2)
%                     disp('ERRROR')
%                     disp(phi)
%                     disp(phi2)
%                     disp(prob{bf})
%                     disp(prob_p)
%                 end
% %                     prob_whole = prob{bf};
%                     prob_p(prob_p==0) = 1; % avoid log0 when computing entropy
%                     H1 = - sum(prob_whole.*log2(prob_p)) ;
% 
%                     prob_whole(prob_whole==0) = 1;
%                     H2 = - sum(prob_whole.*log2(prob_whole));
%                     phi = H1 - H2;
                elseif op_small_phi == 1
                    phi = L1norm(prob{bf},prob_p); 
                elseif op_small_phi == 2
                    phi = emd_hat_gd_metric_mex(prob{bf},prob_p,gen_dist_matrix(length(prob_p)));                  
%                 elseif (op_small_phi == 3)
%                     phi = k_distance(prob{bf},prob_p);
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
             
%             if (all(numerator == [1 2]) && all(denom_past == [1 2])) && (Norm ~= 0)
%                 disp('bf:')
%                 disp(bf)
%                 disp('full dist')
%                 disp(prob{bf})
%                 disp('part dist')
%                 disp(prob_p)
%                 disp('phi:')
%                 disp(phi)
%                 disp('phi_norm')
%                 disp(phi/Norm)
%             end
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
end



function Norm = Normalization(denom_part1,denom_part2,numerator_part1,numerator_part2,xf_1,xf_2)

if nargin == 4
    Norm = min(length(numerator_part1),length(denom_part2)) + min(length(numerator_part2),length(denom_part1));
else %Larissa: this is outdated
    Norm = min(length(numerator_part1),length(denom_part2)) + min(length(numerator_part2),length(denom_part1)) ...
        + min(length(numerator_part1),length(xf_2)) + min(length(numerator_part2),length(xf_1));
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
