function [Big_phi phi_all_values prob_cell MIP M_IRR network] = big_phi_comp_fb(subsystem,whole_sys_state,network)
%%  compute big phi for a subset, subsystem
% subsystem: a subset of the whole system (can be the whole system itself)
% whole_sys_state: current state

% global BRs, global FRs
num_nodes_subsys = length(subsystem);
num_states_subsys = prod([network.nodes(subsystem).num_states]);

op_single = network.options(12);     % just needed for console output
op_console = network.options(8);   %(1: display the results, 0: not)
%% numerator data
% ???This is where we build subsets_subsys of purviews (power-set exclude empty
% set)
subsets_subsys = cell(num_states_subsys-1,1); % subtract one since we don't consider the empty system

% we can do this better for sure - TODO Larissa: Do it as in big_phi_all
% the subsets!!
k = 1;
for subset_size = 1:num_nodes_subsys % can this be done in one for-loop over k = 1:num_states_subsys-1 ?
    C = nchoosek(subsystem,subset_size); % create a matrix of combinations of subsystem of size i
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j = 1:N_C % for all combos of size i
        numerator = C(j,:); % pick a combination
        subsets_subsys{k} = numerator;% store combo
        k = k + 1;
    end
end

% Larissa: Do subsets_subsys the way subsets is done in big_phi_all
% for i = 1:network.num_states-1 % don't include empty set, this is why for-loop starts at 2
%     subsets_subsys{i} = subsystem(logical(network.b_table{i+1,num_nodes_subsys}));
% end
% 

MIP = cell(num_states_subsys-1,1); % MIP in the past, the present, and the future <--?
phi_all_values = zeros(num_states_subsys-1,3); % small phis (for each purview) overall, backward, and forward

prob = cell(num_states_subsys-1,1); % transition repertoire
prob_prod = cell(num_states_subsys-1,1); % partitioned transition repertoire

M_IRR = cell(0,0);

%% computing small phis
EmptyCon = zeros(num_states_subsys-1,1);
for ci=1: num_states_subsys-1  % loop over purview subsets_subsys
    numerator = subsets_subsys{ci}; % given data of numerator
    %Smart purviews: if any element inside the numerator does not have inputs or outputs,
    %no need to calculate purview
    %Larissa: actually one only needs to check WITHIN the Subsystem!!!
    Nconnect = [sum(network.connect_mat(numerator,subsystem),2) sum(network.connect_mat(subsystem,numerator))'];
    EmptyCon(ci) = numel(Nconnect)-nnz(Nconnect);
    %EmptyCon(ci) =0; % Old version without smart purviews, to check
    if EmptyCon(ci) == 0
        [phi_all_values(ci,:) prob{ci} prob_prod{ci} MIP{ci} network] ...
            =  phi_comp_ex(subsystem,numerator,whole_sys_state,subsets_subsys,network);
    else
        phi_all_values(ci,:) = [0 0 0];
        uniform_dist = ones(num_states_subsys,1)/num_states_subsys;
        forward_maxent_dist = comp_pers_cpt(network.nodes,[],subsystem,[],'forward');
        prob{ci} = {uniform_dist, forward_maxent_dist(:)}; 
        prob_prod{ci} = {uniform_dist, forward_maxent_dist(:)}; 
        MIP{ci} = cell(2,2,2); % should we change these to uniform, full sys... etc
    end    
end

prob_cell = cell(2,1);
prob_cell{1} = prob;
prob_cell{2} = prob_prod;
phi = phi_all_values(:,1);
phi_m = zeros(num_nodes_subsys,3); % cumulative sum
index_vec_IRR = find(phi ~= 0);
N_IRR = length(index_vec_IRR); %Number of irreducible concepts

if(N_IRR~=0)
    concepts = zeros(num_states_subsys,N_IRR);  %Larissa: This is only used for strange big phis
    concept_phis = zeros(1,N_IRR);
    j = 1;
    for i = 1:num_states_subsys-1
        if (phi(i) ~= 0)
            concepts(:,j) = prob{i}{1};
            concept_phis(j) = phi(i);
            j = j + 1;
        end
    end

    M_IRR = cell(N_IRR,1);
    for i=1: N_IRR
        j = index_vec_IRR(i);
        M_IRR{i} = subsets_subsys{j};
    end
end

Big_phi = sum(phi); %Now for all big_phi options

% if op_big_phi == 0
%     Big_phi = sum(phi);
% else
%     err = MException('Options:UnSetValue', ...
%         'Option for method of computing Big Phi is incorrect');
%     throw(err);
% end

%% display
if op_console
    for i_C=1: num_states_subsys-1
        C = subsets_subsys{i_C};
        i = length(C);
    
        if abs(phi(i_C)) < 10^-8
            phi(i_C) = 0;
        end
        
        if EmptyCon(i_C) > 0
%             fprintf('C=%s: nodes lack input or output \n',mod_mat2str(C))
        else
            if i > 1 || op_single == 1        
                [string_p string] = make_title_fb(MIP{i_C},0,1); % op_context = 0; op_min = 1
                fprintf('C = %s: Core Concept: %s\n',mod_mat2str(C),string{3});
                fprintf('MIP = %s\n', string_p{3});
                fprintf('\tPast: phi_b = %f\n',phi_all_values(i_C,2));
%                 fprintf('Partition: %s\n',string_p{3});
%                 fprintf('Distribution (full):\n');
%                 disp(prob{i_C}{1}');
%                 fprintf('Distribution (partitioned):\n');
%                 disp(prob_prod{i_C}{1}');
                fprintf('\tFuture: phi_f = %f\n',phi_all_values(i_C,3));
%                 fprintf('Distribution (full):\n');
%                 disp(prob{i_C}{2}');
%                 fprintf('Distribution (partitioned):\n');
%                 disp(prob_prod{i_C}{2}');
            end

            fprintf('\tphi_mip = %f\n\n',phi(i_C));
        end
    end
end
end
