function [phi prob prob_prod_MIP MIP network] = phi_comp_ex(subsystem,numerator,whole_sys_state,subsets_subsys,network)
%% compute small phi for a purview
op_extNodes = network.options(11);

if op_extNodes == 0
    extNodes = setdiff(network.full_system, subsystem);
else
    extNodes = [];
end   

num_nodes_subsys = length(subsystem);
num_states_subsys = prod([network.nodes([subsystem]).num_states]);
   
phi_MIP = zeros(num_states_subsys-1,2);
prob_cand = cell(num_states_subsys-1,1);
prob_prod_MIP_cand = cell(num_states_subsys-1,1);
MIP_cand = cell(num_states_subsys-1,1);

for i=1: num_states_subsys-1
    %Smart purviews: Only test those connections that actually exist
    denom = subsets_subsys{i};
    if nnz(sum(network.connect_mat(numerator,denom),1) == 0) > 0 % some denom is not input of numerator (numerator) --> no phiBR
        if nnz(sum(network.connect_mat(denom,numerator),2) == 0) == 0 % but denom is output
            [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i} network] ...
                = phi_comp_bORf(subsystem,numerator,denom,whole_sys_state,network,2);
        else
            uniform_dist = ones(num_states_subsys,1)/num_states_subsys; % for BR uniform maxent, for FR forward maxent
            forward_maxent_dist = comp_pers_cpt(network.nodes,[],subsystem,whole_sys_state,'forward', extNodes);
            prob_cand{i} = {uniform_dist; forward_maxent_dist(:)};
            prob_prod_MIP_cand{i} = cell(2,1);
            MIP_cand{i} = cell(2,2,2);
        end
    else
        if nnz(sum(network.connect_mat(denom,numerator),2) == 0) > 0 % denom is not output, but denom is input
            [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i} network] ...
                = phi_comp_bORf(subsystem,numerator,denom,whole_sys_state,network,1); 
        else % denom is both
            [phi_MIP(i,:) prob_cand{i} prob_prod_MIP_cand{i} MIP_cand{i} network] ...
                = phi_comp_bf(subsystem,numerator,denom,denom,whole_sys_state,network); 
        end 
    end    
end

%% exlusion principle
max_phi_MIP_bf = zeros(2,1); % backward and forward phi
MIP = cell(2,2,2);
prob = cell(2,1);
prob_prod_MIP = cell(2,1);
for bf = 1:2
    [max_phi_MIP_bf(bf) j_max] = max_ex(phi_MIP(:,bf),subsets_subsys);
    MIP(:,:,bf) = MIP_cand{j_max}(:,:,bf);
    prob{bf} = prob_cand{j_max}{bf};
    prob_prod_MIP{bf} = prob_prod_MIP_cand{j_max}{bf};
    if bf == 1
        xp = subsets_subsys{j_max};
    else
        xf = subsets_subsys{j_max};
    end
end

phi = [0 max_phi_MIP_bf']; % phi = [overall backwards forwards]
phi(1) = min(max_phi_MIP_bf(1),max_phi_MIP_bf(2));

%% imposing maxent on units outside of perspectives
for i = 1:2
    if i == 1
        denom = xp;
        if length(denom) ~= num_nodes_subsys %Larissa: in principle it could happen that although the denominator is smaller prob is already expanded 
                                             %(if there were no connections and it was set to maxent, but
                                             %that should not happen, as it wouldn't j_max
            prob{i} = expand_prob(prob{i},subsystem,denom);
            prob_prod_MIP{i} = expand_prob(prob_prod_MIP{i},subsystem,denom);
        end
    else
        denom = xf;
        if length(denom) ~= num_nodes_subsys 
            denom_rest = pick_rest(subsystem,denom);
            fmaxent_denom_rest = comp_pers_cpt(network.nodes,[],denom_rest,whole_sys_state,'forward',extNodes);
            prob{i} = expand_prob_general(prob{i},subsystem,denom,fmaxent_denom_rest(:));
            prob_prod_MIP{i} = expand_prob_general(prob_prod_MIP{i},subsystem,denom,fmaxent_denom_rest(:));
        end
    end 
end

% if op_console
%     fprintf('Core concept: numerator=%s xp=%s  xf=%s\n',mod_mat2str(numerator),mod_mat2str(xp),mod_mat2str(xf));
%     fprintf('phi=%f\n',phi);
% end
% figure(1)
% subplot(1,2,1),imagesc(prob)
% subplot(1,2,2),imagesc(prob_prod_MIP)
% phi_MIP
% [phi i j] = max2(phi_MIP,subsets_subsys)
% pause;

end

function [X_max i_max j_max] = max2(X,subsets_subsys)
% exclusion principle: if the value is the same, take the bigger one
X_max = -Inf;
i_max = 1;
j_max = 1;
s_max = 0;
for i=1: size(X,1)
    for j=1: size(X,2)
        s = length(subsets_subsys{i}) + length(subsets_subsys{j});
        cond1 = X(i,j) > X_max;
        cond2 = X(i,j) == X_max && s>= s_max;
        if cond1 || cond2
            X_max = X(i,j);
            i_max = i;
            j_max = j;
            s_max = s;
        end
    end
end

end


function [X_max i_max] = max_ex(X,subsets_subsys)
% exclusion principle: if the value is the same, take the bigger one
epsilon = 10^-10;
X_max = -Inf;
i_max = 1;
s_max = 0;
for i=1: size(X,1)
    s = length(subsets_subsys{i});
    cond1 = X(i) > X_max;
    cond2 = abs(X(i) - X_max) < epsilon && s>= s_max;
    if cond1 || cond2
        X_max = X(i);
        i_max = i;
        s_max = s;
    end
end

end