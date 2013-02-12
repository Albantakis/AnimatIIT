function cpt = cpt_factory_tpm(this_node, inputs, nodes, num_total_nodes, tpm)

% THIS FUNCTION CURRENTLY ONLY WORKS FOR BINARY NODES!
% Because it takes the state x nodes tpm

num_sys_nodes = num_total_nodes/2;

dim_sizes = ones(1,num_total_nodes);
dim_sizes(1:num_sys_nodes) = [nodes(1:num_sys_nodes).num_states];

prob_this_node_on = reshape(tpm(:,this_node.num-num_sys_nodes),dim_sizes);


dim_sizes(this_node.num) = this_node.num_states;
cpt = zeros(dim_sizes);
indices = cell(1,num_total_nodes);
indices(:) = {':'};
indices{this_node.num} = 2;

cpt(indices{:}) = prob_this_node_on;

indices = cell(1,num_total_nodes);
indices(:) = {':'};
indices{this_node.num} = 1;

cpt(indices{:}) = 1 - prob_this_node_on;

for i = 1:num_sys_nodes
    
    if ~any(i == inputs)
        cpt = sum(cpt,i)/2;
    end
end
