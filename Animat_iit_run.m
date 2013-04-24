function Animat_iit_run

clear all
%% Animat specificities
numSen = 2;
numMot = 2;
numNodes = 8;
% -------------------------------------------------------------------------
TrialNum =9;
TrialType = 'c35a271_36';
AnimatPath = strcat('~/Documents/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', TrialType , '/trial', int2str(TrialNum),'_');
% -------------------------------------------------------------------------
%% options for one state and KLD and sum of small phis
%in_options = [3     1     2     1     1     0     0     1     1     0     0     0     1     1     1     0     0     1     0];
%all states KLD and sum of small phis
%in_options = [0     1     1     0     0     1     1     0     0     3     1     2     1     0     0     1     1     0];
%all states L1 and L1 norm for complexes
%in_options = [0     1     1     1     1     1     1     0     2     0     0     1     1     0     0     1     1     0];
%all states EMD and L1 norm for constellations only (9: parfor,10:strongconn, 11:freeze)
in_options =  [0     1     0     2     2     1     1     0     1     0     0     1     1     0     0     1     1     0];


   
%Larissa: TODO check this and options in general!
op_average = in_options(2); % 0: use a specific current state 1: average over all possible current states
op_write = 2; % 1 --> write output to file

if op_write > 0
    Foldername = strcat('Freeze_', TrialType,'_trial', int2str(TrialNum));
    mkdir(Foldername)
end 

%% parallel computing
% if a pool is open, close it
if matlabpool('size')
    matlabpool close force;
end
% if parallel option is on, open a new pool
op_parallel = in_options(1);
if op_parallel
    matlabpool;
end

%% begin timer and disp notification
tic
fprintf('\nRunning...\n\n')

%% for loop across generations
for g = [56320       56832       57344];%[30208:512:60000-1 59984];
    results = [];
    
    Animat_gen = g;
        
    J_tempfile = strcat(AnimatPath, int2str(Animat_gen), '_EdgeList.txt');
    J_temp = load(J_tempfile);

    if ~isempty(J_temp)
        J_temp = unique(J_temp, 'rows')+1;
        % MOTORS ARE SET TO 0 -> they don't actually have recurrent
        % connections
        J_temp = J_temp(J_temp(:,1) <= numNodes-numMot,:); 
        %Sensors might have incoming connection, but because they actually don't do anything they shouldn't be taken into account
        J_temp = J_temp(J_temp(:,2) > numSen,:);
        
        J_sparse = sparse(J_temp(:,1), J_temp(:,2),1,numNodes,numNodes);
            
        J = full(J_sparse)';     
        
        %TODO CHECK IF TPM Is correct with large J
        [tpm, used_nodes] = Animat_ReducedTpmSmall(Animat_gen, AnimatPath, numSen, numMot, J);
        tpm(:,used_nodes < numSen+1) = 0.5; %Sensors might be switched on through mechanism but that is overwritten by environment --> doesn't do anything

        
        J = J(used_nodes, used_nodes);
        if ~isempty(J)
            % NODES WITH ONLY OUTPUTS are set to 0 in animat program -> they shouldn't be
            % marginalized -> take them out of TPM
            % STRATEGY: leave motors in TPM but always use 00 state for them
            results.numConn = nnz(J);
            in_connect_mat = J;
            results.connect_mat = J;
            results.used_nodes = used_nodes-1;

            if op_average == 0
                %current_state = zeros(N,1); % all OFF
                %current_state = ones(N,1); % all ON
                current_state = [0.5 0 1 0]';
                state_max = 1;
            else
                [LifeTransitions, p_LifeTransitions] = Animat_LifeStateListSmall(Animat_gen, used_nodes, AnimatPath, numSen);
                PastStates = LifeTransitions(:,1:numel(used_nodes));
                LifeStates = LifeTransitions(:,(numel(used_nodes)+1):2*numel(used_nodes));
                LifeState_index = zeros(1,size(LifeStates,1));
                %state_max = size(LifeStates, 1);
                for i = 1:size(LifeStates,1)
                    LifeState_index(i) = state2index(LifeStates(i,:),2.*ones(length(used_nodes),1));
                end    
                results.LifeState_index = LifeState_index;    
                results.p_LifeTransitions = p_LifeTransitions;
            end

    %% new way of input with nodes
            num_nodes = size(tpm,2);
            logic_types = cell(num_nodes,1);

            % init struct array
            %Larissa: Logic_type has to be empty for the removals to work
            %correctly!!!
            in_nodes(2*num_nodes) = struct('num',2*num_nodes,'name',[num2str(used_nodes(num_nodes)-1) '_c'],'num_states',2,...
                                        'state_names',{{'0' '1'}},'logic_type',logic_types{num_nodes},'cpt',[],...
                                        'num_sys_nodes',num_nodes,'input_nodes',[]);

            % make past node structs                        
            for i = 1:num_nodes

                in_nodes(i) = struct('num',i,'name',[num2str(used_nodes(i)-1) '_p'],'num_states',2,...
                                        'state_names',{{'0' '1'}},'logic_type',logic_types{i},'cpt',[],...
                                        'num_sys_nodes',num_nodes,'input_nodes',[]);

            end

            % make current node structs and their tpms
            for i = 1:num_nodes

                in_nodes(num_nodes + i) = struct('num',num_nodes + i,'name',[num2str(used_nodes(i)-1) '_c'],'num_states',2,...
                                        'state_names',{{'0' '1'}},'logic_type',logic_types{i},'cpt',[],...
                                        'num_sys_nodes',num_nodes,'input_nodes',[]);

                input_nodes = 1:num_nodes;
                input_nodes_indices = input_nodes(logical(in_connect_mat(i,:)));
                in_nodes(num_nodes + i).input_nodes = input_nodes_indices;

            %     input_nodes = nodes(input_nodes_indices);
            %     nodes(num_nodes + i).cpt = cpt_factory_mechs(nodes(num_nodes + i),input_nodes,2*num_nodes,noise);
            %     disp(nodes(num_nodes + i).cpt)
            %     test_cpt = cpt_factor_tpm(nodes(num_nodes + i), input_nodes_indices, nodes, 2*num_nodes, tpm);
                in_nodes(num_nodes + i).cpt = cpt_factory_tpm(in_nodes(num_nodes + i), input_nodes_indices, in_nodes, 2*num_nodes, tpm);
            end

    %% initialize data which is the same for all states

            network.connect_mat = in_connect_mat;
            network.options = in_options;
            network.nodes = in_nodes;
            network.num_nodes = num_nodes;
            network.tpm = tpm;
            network.full_system = 1:num_nodes;
            network.num_subsets = 2^num_nodes;
            network.noise = 0;
            network.num_states = prod([network.nodes(network.full_system).num_states]);
            % Larissa: check, cause we should use lifestates, but maybe all network
            % states is actually needed some time later
%             network.states = zeros(network.num_nodes,network.num_states);
%             for i = 0:network.num_states - 1
%                 network.states(:,i+1) = dec2multibase(i,[network.nodes(network.full_system).num_states]);
%             end

            % Larissa: I don't think this changes, so can be put here for all
            % states
            % binary table and states list
            % need to rethink use of b_table when allowing for more than binary nodes
            network.b_table = cell(network.num_subsets,network.num_nodes);
            for i = network.full_system
                for j = 1:2^i
                    network.b_table{j,i} = trans2(j-1,i);
                end
            end

            if op_average ~= 0
                %state_max = network.num_states;
                state_max = numel(LifeState_index);
            end    
            % init cell arrays for results - OLD WAY
            Big_phi_M_st = cell(state_max,1);
            Big_phi_MIP_st = cell(state_max,1);
            MIP_st = cell(state_max,1);
            Complex_st = cell(state_max,1);
            BFCut_st = cell(state_max,1);
            prob_M_st = cell(state_max,1);
            phi_M_st = cell(state_max,1);
            concept_MIP_M_st = cell(state_max,1);
            complex_MIP_M_st = cell(state_max,1);
            Big_phi_MIP_all_M_st = cell(state_max,1);
            complex_MIP_all_M_st = cell(state_max,1);
            purviews_M_st = cell(state_max,1);

            % find main complex (do system partitions)
            op_complex = network.options(3);
            %% main loop over life states
            for z = 1:state_max
                if op_average == 1
                    %current_state = network.states(:,z);
                    current_state = LifeStates(z, :)';
                    past_state = PastStates(z,:)';
                end
                this_state = current_state;
                network.current_state = current_state;
                network.past_state = past_state;
                % init backward rep and forward reps for each state
                network.BRs = cell(network.num_subsets); % backward repertoire
                network.FRs = cell(network.num_subsets); % forward repertoire

                %fprintf(['State: ' num2str(this_state') '\n'])

                %Larissa: change that to only calculate Lifestates
                % is it possible to reach this state
                check_prob = partial_prob_comp(network.full_system,network.full_system,this_state,tpm,network.b_table,1); % last argument is op_fb = 1;
                state_reachable = sum(check_prob);

                if ~state_reachable %|| ~ismember(z, LifeState_index)   %Shouldn't happen!

                    %fprintf('\tThis state cannot be realized...\n')

                    Big_phi_M_st{z} = NaN;
                    Big_phi_MIP_st{z} = NaN;

                    % SET OTHERS

                else
                    % only consider whole system
                    % THIS OPTION NEEDS TO BE WORKED OUT!
                    if op_complex == 0 

            %             [BRs FRs] = comp_pers(this_state,tpm,b_table,options);
            %             (network.full_system,this_state,tpm,network.b_table,network.options)
                        % [Big_phi phi prob_cell MIPs M_IRR network] = big_phi_comp_fb(network.full_system,this_state,network);
                        % Big_phi_M_st{z} = Big_phi;
                        [Big_phi_M phi_M prob_cell concept_MIP_M purviews_M] = big_phi_comp_fb(network.full_system,this_state,network);
                        
                        Complex_st{z} = network.full_system;
                        
                        if op_write == 2
                            results.state(z) = Animat_rewrap_ZombieData(network.full_system, Big_phi_M, phi_M, concept_MIP_M, purviews_M);
                        end        
                    % find the complex
                    elseif op_complex == 1

                        [Big_phi_M phi_M prob_M M_cell concept_MIP_M purviews_M network concept_MIP_M_subs] = big_phi_all(network, this_state);                                             
                        % complex search
                        [Big_phi_MIP MIP Complex M_i_max BFCut Big_phi_MIP_M complex_MIP_M Big_phi_MIP_all_M complex_MIP_M_all BFCut_M] = ...
                            complex_search(Big_phi_M,M_cell, purviews_M, network.num_nodes,prob_M,phi_M,network.options,concept_MIP_M,network);
    
        
                        Big_phi_M_st{z} = Big_phi_M;
                        Big_phi_MIP_st{z} = Big_phi_MIP_M;
                        % it looks like MIP is never  used
                        MIP_st{z} = MIP;
                        Complex_st{z} = Complex;
                        prob_M_st{z} = prob_M;
                        phi_M_st{z} = phi_M;
                        BFCut_st{z} = BFCut; %M1->M2 noised, or M1<-M2
                        
                        concept_MIP_M_st{z} = concept_MIP_M;
                        complex_MIP_M_st{z} = complex_MIP_M;
                        Big_phi_MIP_all_M_st{z} = Big_phi_MIP_all_M;
                        complex_MIP_all_M_st{z} = complex_MIP_M_all;
                        purviews_M_st{z} = purviews_M;

                        if op_write == 1
                            results.state(z) = rewrap_data(Big_phi_M, phi_M, prob_M, M_cell, concept_MIP_M, purviews_M,...
                                    Big_phi_MIP, MIP, Complex, M_i_max,  Big_phi_MIP_M, complex_MIP_M, Big_phi_MIP_all_M, complex_MIP_M_all);
                        elseif op_write == 2
                            results.state(z) = Animat_rewrap_data_short(Big_phi_M, phi_M, prob_M, M_cell, concept_MIP_M, purviews_M,...
                                    Big_phi_MIP, MIP, Complex, M_i_max, BFCut);
                        end    
                    end
                end
            end        

        %% store output data

        %Larissa: Make one file for all states that contains everything except the
        %largest part, for the largest part make one output per state
        if op_write == 1
            results.network = network;

            output_data.tpm = tpm;
            output_data.J = network.connect_mat;
            output_data.current_state = current_state;
            output_data.noise = 0;
            output_data.options = network.options;
            output_data.num_nodes = num_nodes;
            output_data.Big_phi_M = Big_phi_M_st;
            output_data.Big_phi_MIP = Big_phi_MIP_st;
            % KILL THIS ONE BELOW
            output_data.MIP = MIP_st;
            output_data.Complex = Complex_st;
            output_data.concepts_M = prob_M_st;
            output_data.small_phi_M = phi_M_st;
            output_data.concept_MIP_M = concept_MIP_M_st;
            output_data.complex_MIP_M = complex_MIP_M_st;
            output_data.M_cell = M_cell;
            output_data.Big_phi_MIP_all_M = Big_phi_MIP_all_M_st;
            output_data.complex_MIP_all_M = complex_MIP_all_M_st;
            output_data.purviews_M = purviews_M_st;
        end

        results.Complex = Complex_st;
        results.options = in_options;

        %% finish & cleanup: stop timer, save data, open explorer gui, close matlabpool
        toc

        if op_write == 1
           cd(Foldername)
           StateFilename = strcat('Animat_gen', int2str(Animat_gen),'_Results');      
           save(StateFilename, 'results');
           Filename = strcat('Animat_gen', int2str(Animat_gen),'_GuiInput');      
           save(Filename, 'output_data');
           cd ..
        elseif op_write == 2
           cd(Foldername)
           if op_complex == 0
               StateFilename = strcat('Animat_gen', int2str(Animat_gen),'_Zombie');
           elseif op_complex == 1
               StateFilename = strcat('Animat_gen', int2str(Animat_gen),'_ShortResults');
           end
           save(StateFilename, 'results');
           cd ..
        end  

        %save('last_run_output.mat','output_data','-v7.3');
        % save('save_test1.mat','output_data');
        % save('save_test2.mat','output_data','-v7.3');

        end
    end
end % for loop g, generations

if op_parallel
    matlabpool close force;
end


