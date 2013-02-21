load GT_standardExample
% iit_run(tpm, in_connect_mat, current_state, in_noise, in_options, in_nodes)
connectivity_matrix = connectivity_mat;
clear nodes;

% op_parallel = in_options(1);
% op_average = network.options(2);
% op_complex = network.options(3);
% op_small_phi = network.options(4); % 0 KLD 1 L1 2 EMD
% op_big_phi = network.options(5); % 0 KLD 1 L1 2 EMD
% op_normalize = network.options(6); % small phi
% op_normalize = options(7);         % big phi
% op_console = options(8);
% op_single = network.options(11);    % just needed for console output == 1


%one state KLD and sum of small phis 
in_options = [0     0     1     0     0     1     1     0  0 3     1     2     1     0     0     1     1     0]; 


%one state L1 and L1 norm for complexes
% in_options = [0     1     1     1     1     1     1     0     0  3     1     2     1     0     0     1     1     0];


num_nodes = 3;
nodes(2*num_nodes) = struct('num',2*num_nodes,'name',[num2str(num_nodes) '_c'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types(num_nodes),'cpt',[],...
                            'num_sys_nodes',num_nodes,'input_nodes',[]);
for i = 1:num_nodes
    
    nodes(i) = struct('num',i,'name',[num2str(i) '_p'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types(i),'cpt',[],...
                            'num_sys_nodes',num_nodes,'input_nodes',[]);
    
end

% make current node structs and their tpms
for i = 1:num_nodes
    
    nodes(num_nodes + i) = struct('num',num_nodes + i,'name',[num2str(i) '_c'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types(i),'cpt',[],...
                            'num_sys_nodes',num_nodes,'input_nodes',[]);

	input_nodes = 1:num_nodes;
    input_nodes_indices = input_nodes(logical(connectivity_matrix(i,:)));
    nodes(num_nodes + i).input_nodes = input_nodes_indices;
    
%     input_nodes = nodes(input_nodes_indices);
%     nodes(num_nodes + i).cpt = cpt_factory_mechs(nodes(num_nodes + i),input_nodes,2*num_nodes,noise);
%     disp(nodes(num_nodes + i).cpt)
%     test_cpt = cpt_factory_tpm(nodes(num_nodes + i), input_nodes_indices, nodes, 2*num_nodes, tpm);
    nodes(num_nodes + i).cpt = cpt_factory_tpm(nodes(num_nodes + i), input_nodes_indices, nodes, 2*num_nodes, tpm);
    
%     if any(nodes(num_nodes + i).cpt ~= test_cpt)
%         disp('error')
%     end
    
    
    
end

iit_run(tpm, connectivity_matrix, cur_state', 0, in_options, nodes)
