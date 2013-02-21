function Rem_network = cpt_removal_network(this_subset, network)
% build a cell array that contains all of the subsets
N = network.num_nodes;
Nsub = numel(this_subset);
prev_logic_gates = zeros(Nsub,1);
prev_inputs = cell(Nsub,1);
new_input = cell(Nsub,1);
inputs_sub = cell(Nsub,1);

subsys_nodes = 1:Nsub;

Jsub = network.connect_mat(this_subset, this_subset);
for i = 1:Nsub
        inputs_sub{i} = subsys_nodes(logical(Jsub(i,:)));
end
    
if isempty(network.nodes(1).logic_type) %if tpm is loaded and logic_types are not defined
    new_logic_gates = cell(Nsub,1); 
 
    states = network.states';
    
    %Exclude Only Output Nodes
    excluded_nodes = setdiff(1:N,this_subset);
    %take TPM for all excluded elements = 0. 
    ind0 = find(sum(states(:,excluded_nodes),2) == 0);
    new_tpm = network.tpm(ind0,this_subset);   
else    
    % Check remaining inputs to each node and redefine logic gates
    for i = 1:numel(this_subset)
        prev_logic_gates(i) = network.nodes(this_subset(i)).logic_type;
        prev_inputs{i} = network.nodes(this_subset(i)+network.num_nodes).input_nodes;
        new_input{i} = this_subset(inputs_sub{i});
    end

    %Larissa -> Should be 0 now is es aber nicht!
    new_logic_gates = num2cell(convert_logic_gates_removal(prev_logic_gates, prev_inputs, new_input)); %this becomes a cell to be able to make it [] in the case without logic_gates
    
    new_tpm = zeros(2^Nsub,Nsub);
    % new TPM
    for k = 1:2^Nsub   
        x0 = trans2(k-1,Nsub);
        for i = 1:Nsub
            input_vec = x0(inputs_sub{i});  
            new_tpm(k,i) = logic_gates(input_vec,new_logic_gates{i},network.noise);
        end
    end
end

% setup node strucs and cpts for each node
% inputs = struct('num',{1 2},'name',{'A_p' 'B_p'},'num_states',{2 2},'state_names',{{'0' '1'}},'logic_type',{2 3})
logic_types = new_logic_gates;
% init struct array
nodes(2*Nsub) = struct('num',2*Nsub,'name',[num2str(Nsub) '_c'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types{Nsub},'cpt',[],...
                            'num_sys_nodes',Nsub,'input_nodes',[]);
% make past node structs                        
% Larissa: Not sure yet if for num it should just be 'i' or 'this_subset(i)'
for i = 1:Nsub
    nodes(i) = struct('num',i,'name',[num2str(this_subset(i)) '_p'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types{i},'cpt',[],...
                            'num_sys_nodes',Nsub,'input_nodes',[]);
end

% make current node structs and their tpms
for i = 1:Nsub
    nodes(Nsub + i) = struct('num',Nsub + i,'name',[num2str(this_subset(i)) '_c'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types{i},'cpt',[],...
                            'num_sys_nodes',Nsub,'input_nodes',[]);
    nodes(Nsub + i).input_nodes = inputs_sub{i};

    nodes(Nsub + i).cpt = cpt_factory_tpm(nodes(Nsub + i), inputs_sub{i}, nodes, 2*Nsub, new_tpm);
end

Rem_network.this_subset = this_subset;
Rem_network.connect_mat = Jsub;
Rem_network.options = network.options;
Rem_network.nodes = nodes;
Rem_network.num_nodes = Nsub;
Rem_network.tpm = new_tpm;
Rem_network.full_system = subsys_nodes;
Rem_network.num_subsets = 2^Nsub;
Rem_network.current_state = network.current_state(this_subset);
Rem_network.past_state = [];
Rem_network.num_states = prod([Rem_network.nodes(subsys_nodes).num_states]);
Rem_network.noise = network.noise;
Rem_network.b_table = network.b_table(1:2^Nsub,1:Nsub); %just used in phi_comp_ex_unidir, better replace by set of subsets
% Rem_network.states =      % probably never used
Rem_network.BRs = cell(2^Nsub);
Rem_network.FRs = cell(2^Nsub);
end
