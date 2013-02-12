function iit_run(tpm, in_J, current_state, in_noise, in_options, in_nodes)
% IIT_RUN Computes concept, small and big phi, and partition information
% for all subsets of a system (exluding the empty set) over a binary network.
%
%   IIT_RUN(TPM, J, CURRENT_STATE, NOISE, OPTIONS) takes a TPM in
%   state x node form, that is TPM(i,j) is the probability node_i = 1 at 
%   time t+1 given that the system is in state j at time t. J is the
%   connectivity matrix of the network such that J(i,j) = 1 when j has a
%   directed edge to i, and J(i,j) = 0 otherwise. current_state is the
%   state of the system at time t (only used if the options are not set to
%   compute over all states). NOISE is a global noise put on all
%   outgoing messages and must take a value on the interval [0,.5]. OPTIONS
%   is a structure of options for the algoirthm created using the
%   set_options function
%
%   see also set_options

% parallel computing
if matlabpool('size')
    matlabpool close force;
end

op_parallel = in_options(19);

if op_parallel
    matlabpool;
end

% profile on

tic

fprintf('\nRunning...\n\n')

% get num_nodes, the number of nodes in the whole system
num_nodes = length(in_nodes)/2;

% global inputs
 
% global 

% global BRs_check FRs_check
% global BRs_check2 FRs_check2

% global func_time inline_time cpt_time tpm_time
% func_time = 0;
% inline_time = 0;
% cpt_time = 0;
% tpm_time = 0;

network.J = in_J;
network.options = in_options;
network.nodes = in_nodes;
network.num_nodes = num_nodes;
network.tpm = tpm;
% options(10) = 1;

output_data.tpm = tpm;
output_data.J = network.J;
output_data.current_state = current_state;
output_data.noise = in_noise;
output_data.options = network.options;
output_data.num_nodes = num_nodes;

network.full_system = 1:num_nodes;
network.num_subsets = 2^num_nodes;
network.num_states = prod([network.nodes(network.full_system).num_states]);

% binary table and states list
% need to rethink use of b_table when allowing for more than binary nodes
network.b_table = cell(network.num_subsets,network.num_nodes);
network.states = zeros(network.num_nodes,network.num_states);
for i = network.full_system
    for j = 1:2^i
        network.b_table{j,i} = trans2(j-1,i); % CONSIDER FLIPPING THIS LR
%         if i == network.num_nodes
%             network.states(:,j) = trans2(j-1,i);
%         end
    end
end

for i = 0:network.num_states - 1
    network.states(:,i+1) = dec2multibase(i,[network.nodes(network.full_system).num_states]);
end

% determine if we are computing over all states or just one
op_ave = network.options(18);
if op_ave == 0
    state_max = 1;
    network.states(:,1) = current_state;
else
    state_max = network.num_subsets;
end

% we should deal with different arguments not being included
% if nargin == 4 
%     J = ones(num_nodes);
% elseif nargin == 5
%     J = in_J;
% end


output_data.states = network.states;


% find main complex (do system partitions)
op_complex = network.options(15);

% init cell arrays for results
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


% for each state
for z = 1:state_max
    
    this_state = network.states(:,z);
    
    % init backward rep and forward reps for each state
    network.BRs = cell(network.num_subsets); % backward repertoire
    network.FRs = cell(network.num_subsets); % forward repertoire
    
%     [BRs_check2 FRs_check2] = comp_pers(this_state,tpm,b_table,options);
    
    fprintf(['State: ' num2str(this_state') '\n'])
   
    % is it possible to reach this state
    check_prob = partial_prob_comp(network.full_system,network.full_system,this_state,tpm,network.b_table,1); % last argument is op_fb = 1;
    state_check1 = sum(check_prob);
    
%     tic
%     state_check2 = all(tpm(logical(this_state),:) > 0) & all(tpm(pick_rest(full_system,full_system(this_state)),:) < 1)
%     toc
    
%     disp (state_check1 == state_check2)

    if state_check1 == 0
        
        fprintf('\tThis state cannot be realized...\n')
        
        Big_phi_M_st{z} = NaN;
        Big_phi_MIP_st{z} = NaN;
        
        % SET OTHERS
        
    else
        
        fprintf('\tComputing state...\n')
        
        % only consider whole system
        % THIS OPTION NEEDS TO BE WORKED OUT!
        if op_complex == 0 

%             [BRs FRs] = comp_pers(this_state,tpm,b_table,options);
%             (network.full_system,this_state,tpm,network.b_table,network.options)
            [Big_phi phi prob_cell MIPs M_IRR network] = big_phi_comp_fb(network.full_system,this_state,network);
            % irreducible points
%             [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_cell,phi,MIPs,M, 0,op_fb);
%             plot_REP(Big_phi, IRR_REP,IRR_phi,IRR_MIP, 1, M, options)


            Big_phi_M_st{z} = Big_phi;
            
            % TODO: WE NEED TO HANDLE MIP IN THIS CASE EVEN WE DON'T FIND
            % THE COMPLEX
            
        % find the complex    
        elseif op_complex == 1 
            
            
%             [MIP Complex Big_phi_M Big_phi_MIP_M prob_M phi_M...
%                 concept_MIP_M complex_MIP_M M_cell Big_phi_MIP_all_M complex_MIP_M_all purviews_M] ...
%                 = big_phi_complex(this_state,tpm);
%             
%             
%             [MIP Complex Big_phi_M Big_phi_MIP_M prob_M phi_M concept_MIP_M complex_MIP_M M_cell Big_phi_MIP_all_M complex_MIP_M_all M_IRR_M]
%             big_phi_all(this_state,tpm,options)
            [Big_phi_M phi_M prob_M M_cell concept_MIP_M purviews_M] = big_phi_all(network, this_state);
                                                                
            % complex search
            [Big_phi_MIP MIP Complex M_i_max  Big_phi_MIP_M complex_MIP_M Big_phi_MIP_all_M complex_MIP_M_all] = ...
                complex_search(Big_phi_M,M_cell, purviews_M, network.num_nodes,prob_M,phi_M,network.options);
            
            
            assignin('base','Big_Phi_M',Big_phi_M)
            assignin('base','phi_M',phi_M)
            assignin('base','prob_M',prob_M)
            assignin('base','M_cell',M_cell)
            assignin('base','concept_MIP_M',concept_MIP_M)
            assignin('base','purviews_M',purviews_M)
            assignin('base','Big_phi_MIP',Big_phi_MIP)
            assignin('base','MIP',MIP)
            assignin('base','Complex',Complex)
            assignin('base','M_i_max',M_i_max)
            assignin('base','Big_phi_MIP_M',Big_phi_MIP_M)
            assignin('base','complex_MIP_M',complex_MIP_M)
            assignin('base','Big_phi_MIP_all_M',Big_phi_MIP_all_M)
            assignin('base','complex_MIP_M_all',complex_MIP_M_all)
            
            
            
            Big_phi_M_st{z} = Big_phi_M;
            Big_phi_MIP_st{z} = Big_phi_MIP_M;
            MIP_st{z} = MIP;
            Complex_st{z} = Complex;
            prob_M_st{z} = prob_M;
            phi_M_st{z} = phi_M;
            concept_MIP_M_st{z} = concept_MIP_M;
            complex_MIP_M_st{z} = complex_MIP_M;
            Big_phi_MIP_all_M_st{z} = Big_phi_MIP_all_M;
            complex_MIP_all_M_st{z} = complex_MIP_M_all;
            purviews_M_st{z} = purviews_M;
                

        end
    end

end

% load up output_data struct
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


% if op_ave == 1
%     if op_fb == 0
%         Big_phi_ave = sum(Big_phi_st)/2^num_nodes;
%     else
%         Big_phi_ave = sum(p_x1 .* Big_phi_st); %weighted ave/expected value
%     end
% %     fprintf('Big_phi_ave=%f\n',Big_phi_ave);
% end
toc


% profile off
% profile viewer

fprintf('Loading GUI... \n');

save('output_data_sample.mat','output_data');

iit_explorer(output_data)

if op_parallel
    matlabpool close force;
end


% disp('FUNCTION TIME:')
% disp(func_time)
% disp('INLINE TIME:')
% disp(inline_time)
% disp('CPT TIME:')
% disp(cpt_time)
% disp('TPM TIME:')
% disp(tpm_time)




% mpiprofile viewer