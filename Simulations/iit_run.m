function iit_run(tpm, in_connect_mat, current_state, in_noise, in_options, in_nodes, past_state)
% IIT_RUN Computes concepts, small and big phi, and partition information
% for all subsets of a system (exluding the empty set) over a binary network.
%
%   IIT_RUN(TPM, connect_mat, CURRENT_STATE, NOISE, OPTIONS) takes a TPM in
%   state X node form, that is TPM(i,j) is the probability node_i = 1 at 
%   time t+1 given that the system is in state j at time t. connect_mat is the
%   connectivity matrix of the network such that connect_mat(i,j) = 1 when j has a
%   directed edge to i, and connect_mat(i,j) = 0 otherwise. current_state is the
%   state of the system at time t (only used if the options are not set to
%   compute over all states). NOISE is a global noise put on all
%   outgoing messages and must take a value on the interval [0,.5]. OPTIONS
%   is a structure of options for the algoirthm created using the
%   set_options function
%
%   see also set_options

if nargin < 7
    past_state = [];
end    
%% parallel computing
% in_options(9) = 0;
% if a pool is open, close it
if matlabpool('size')
    matlabpool close force;
end

% if parallel option is on, open a new pool
op_parallel = in_options(1);
op_PHIconcept_fig = 0;
op_extNodes = in_options(11);
if op_parallel
    matlabpool;
end
%% begin timer and disp notification
tic

fprintf('\nRunning...\n\n')

%% initialize data

% get num_nodes, the number of nodes in the whole system
% note that in_nodes is the number of nodes in the GRAPH = 2*num_nodes
num_nodes = length(in_nodes)/2;

network.connect_mat = in_connect_mat;
network.options = in_options;
network.nodes = in_nodes;
network.num_nodes = num_nodes;
network.tpm = tpm;
network.full_system = 1:num_nodes;
network.num_subsets = 2^num_nodes;
network.current_state = current_state;
network.past_state = past_state;
network.num_states = prod([network.nodes(network.full_system).num_states]);

% get rid of everyting below
output_data.tpm = tpm;
output_data.J = network.connect_mat;
output_data.current_state = current_state;
output_data.noise = in_noise;
output_data.options = network.options;
output_data.num_nodes = num_nodes;

% output_data.tpm = tpm;
% output_data.current_state = current_state;
network.noise = in_noise;
% output_data.num_nodes = num_nodes;

% binary table and states list
% need to rethink use of b_table when allowing for more than binary nodes
network.b_table = cell(network.num_subsets,network.num_nodes);
for i = network.full_system
    for j = 1:2^i
        network.b_table{j,i} = trans2(j-1,i);
    end
end

network.states = zeros(network.num_nodes,network.num_states);
for i = 0:network.num_states - 1
    network.states(:,i+1) = dec2multibase(i,[network.nodes(network.full_system).num_states]);
end

%% setup main function call
% determine if we are averaging over all states or just one
op_average = network.options(2);
if op_average == 0
    state_max = 1;
%     network.states(:,1) = current_state;
else
    state_max = network.num_states;
end

% we should deal with different arguments not being included
% if nargin == 4 
%     connect_mat = ones(num_nodes);
% elseif nargin == 5
%     connect_mat = in_connect_mat;
% end

% find main complex (do system partitions)
op_complex = network.options(3);


% init output structs - NEW WAY!
% output_data.results.state(state_max).subsystem.Phi = 0;


% init cell arrays for results - OLD WAY
Big_phi_M_st = cell(state_max,1);
Big_phi_MIP_st = cell(state_max,1);
MIP_st = cell(state_max,1);
Complex_st = cell(state_max,1);
prob_M_st = cell(state_max,1);
phi_M_st = cell(state_max,1);
concept_MIP_M_st = cell(state_max,1);
complex_MIP_M_st = cell(state_max,1);
Big_phi_MIP_all_M_st = cell(state_max,1);
complex_MIP_all_M_st = cell(state_max,1);
purviews_M_st = cell(state_max,1);
BFCut_st = cell(state_max,1);
BFCut_M_st = cell(state_max,1);
M_cell = cell(network.num_subsets-1,1);

%% main loop over states

% for each state
for z = 1:state_max
      
    if op_average
        this_state = network.states(:,z);
    else
        this_state = current_state;
    end
    
    % init backward rep and forward reps for each state
    network.BRs = cell(network.num_subsets); % backward repertoire
    network.FRs = cell(network.num_subsets); % forward repertoire
        
    fprintf(['State: ' num2str(this_state') '\n'])
   
    % is it possible to reach this state
    check_prob = partial_prob_comp(network.full_system,network.full_system,this_state,tpm,network.b_table,1); % last argument is op_fb = 1;
    state_reachable = sum(check_prob);
    
    if ~state_reachable
        
        fprintf('\tThis state cannot be realized...\n')
        
        Big_phi_M_st{z} = NaN;
        Big_phi_MIP_st{z} = NaN;
    
    else
        
        fprintf('\tComputing state...\n')      
 
        if op_complex == 0 %Larissa: Quick and dirty fix, so that it can be loaded into GUI

            % compute big phi in every possible subset
            Big_phi_M = zeros(network.num_states-1,1); % Big_phi for each subset except the empty set
            phi_M = cell(network.num_states-1,1);
            prob_M = cell(network.num_states-1,2); 
            concept_MIP_M = cell(network.num_states-1,1); % the partition that gives Big_phi_MIP for each subset
            purviews_M = cell(network.num_states-1,1);
            M_cell= cell(network.num_states-1,1);
            
            M_cell{end} = network.full_system;
            
            [Big_phi_M(end) phi_M{end} prob_cell concept_MIP_M{end} purviews_M{end}] = big_phi_comp_fb(network.full_system,this_state,network);
            
            % concept distributions
            prob_M(end,:) = prob_cell(:); % first layer is subset, second is purview, third is backward/forward  
            
        % find the complex
        elseif op_complex == 1
           
            [Big_phi_M phi_M prob_M M_cell concept_MIP_M purviews_M network concept_MIP_M_subs] = big_phi_all(network, this_state); %Larissa: this_state should be obsolete as it is in network
            
        end                                                      
        % complex search
        [Big_phi_MIP MIP Complex M_i_max BFCut Big_phi_MIP_M complex_MIP_M Big_phi_MIP_all_M complex_MIP_M_all BFCut_M] = ...
            complex_search(Big_phi_M,M_cell, purviews_M, network.num_nodes,prob_M,phi_M,network.options,concept_MIP_M,network);

        Big_phi_M_st{z} = Big_phi_M;
%             output_data.results.state(z).Phi = Big_phi_M; 
        Big_phi_MIP_st{z} = Big_phi_MIP_M;
%             output_data.results.state(z).Phi_MIP = Phi_MIP;
        % it looks like MIP is never used
        MIP_st{z} = MIP;

        Complex_st{z} = Complex;
%             output_data.results(z).complex_set = complex_set;

        prob_M_st{z} = prob_M;
%             output_data.results(z).concepts = concepts;

        phi_M_st{z} = phi_M;

        BFCut_st{z} = BFCut; %M1->M2 noised, or M1<-M2
        BFCut_M_st{z} = BFCut_M;

        % For removals, the concepts don't yet have the right node
        % names
        if op_extNodes == 1
           for i = 1:size(Big_phi_M,1)-1 %all except full system
               if ~isempty(network.removal_networks{i})
                    this_subset = network.removal_networks{i}.this_subset;
                    for j = 1:size(purviews_M{i},1)
                        purviews_M{i}{j} = this_subset(purviews_M{i}{j});
                    end    
                end
           end
           concept_MIP_M = {concept_MIP_M_subs{1:end-1} concept_MIP_M{end}};
         end    

        concept_MIP_M_st{z} = concept_MIP_M;
        complex_MIP_M_st{z} = complex_MIP_M;
        Big_phi_MIP_all_M_st{z} = Big_phi_MIP_all_M;
        complex_MIP_all_M_st{z} = complex_MIP_M_all;
        purviews_M_st{z} = purviews_M;

        %BFcut not in rewrap_data, but then we need to restructure this
        %anyways
        %output_data.results.state(z) = rewrap_data(Big_phi_M, phi_M, prob_M, M_cell, concept_MIP_M, purviews_M,...
        %            Big_phi_MIP, MIP, Complex, M_i_max,  Big_phi_MIP_M, complex_MIP_M, Big_phi_MIP_all_M, complex_MIP_M_all);

        if op_PHIconcept_fig ==1 
            [CutDistr] = PHI_Cut_concepts(Complex,MIP{1},BFCut,purviews_M, prob_M, phi_M,concept_MIP_M, network); 
        end                           
    end
end

%% store output data

output_data.network = network;

output_data.Big_phi_M = Big_phi_M_st;
output_data.Big_phi_MIP = Big_phi_MIP_st;
output_data.BFCut = BFCut_st;
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
output_data.BFCut_M = BFCut_M_st;


%% finish & cleanup: stop timer, save data, open explorer gui, close matlabpool
toc

fprintf('Loading GUI... \n');

%The tag is only necessary for large networks and then it is very big anyways
%save('last_run_output.mat','output_data','-v7.3'); 
save('last_run_output.mat','output_data');
% save('save_test1.mat','output_data');
% save('save_test2.mat','output_data','-v7.3');

iit_explorer(output_data)

if op_parallel
    matlabpool close force;
end