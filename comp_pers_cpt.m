function perspective = comp_pers_cpt(nodes,num_nodes_indices,denom_nodes_indices,numerator_state,bf_option,extNodes,past_state, M1, M2, bfcut_option)

%  compute BRs and FRs for a single perspective but given some fixed
%  current state

if nargin < 10
    M1 = []; M2 = []; bfcut_option = [];
end
if nargin < 7
    past_state = [];
end    
if nargin < 6
    extNodes = [];
end    

if isempty(denom_nodes_indices)
    perspective = [];
    return
% elseif isempty(num_nodes_indices)
% %     num_sys_nodes = denom_nodes_indices(1).num_sys_nodes;
% %     denom_conditional_joint_size = ones(1,2*num_sys_nodes);
% %     denom_conditional_joint_size(1:num_sys_nodes == denom_nodes_indices
%     denom_conditional_joint = [];
%     return
end

num_sys_nodes = nodes(1).num_sys_nodes;

if strcmp(bf_option,'backward')
    
    denom_nodes = nodes(denom_nodes_indices);
    num_nodes_shift = num_nodes_indices + num_sys_nodes;
    numerator_nodes = nodes(num_nodes_shift);
    
    % no nodes in numerator means maxent over denom
    if isempty(num_nodes_indices)
        
        perspective_dim_sizes = ones(1,num_sys_nodes);
        perspective_dim_sizes(denom_nodes_indices) = [denom_nodes.num_states];
        perspective = ones([perspective_dim_sizes, 1])./prod(perspective_dim_sizes);    %The additional 1 is to take care of selfloops.
        return
        
    end
    
    % this just defines the final dimension of the distribution
    numerator_conditional_joint_size = ones(1,2*num_sys_nodes);
    numerator_conditional_joint_size(denom_nodes_indices) = [denom_nodes.num_states];
    numerator_conditional_joint = ones(numerator_conditional_joint_size);
    
    % setup cell array for conditioning
    conditioning_indices = cell(1,2*num_sys_nodes);
    conditioning_indices(:) = {':'};

    prob_current_state = 1;
    
    % Choose dimensions that are congruent with current state (dim 1 if OFF, dim 2 if ON)
    for i = 1:length(num_nodes_indices) %Loop over numerator nodes
        
        this_node_conditioning_indices = conditioning_indices;
        this_node_conditioning_indices{numerator_nodes(i).num} = numerator_state(numerator_nodes(i).num - num_sys_nodes) + 1;
        next_num_node_distribution = numerator_nodes(i).cpt(this_node_conditioning_indices{:});
        %Larissa: This doesn't seem to do anything...
        prob_current_state = prob_current_state * sum(next_num_node_distribution(:));

        % marginalize over nodes not in denom, these nodes are outside the
        % system for this iteration or they are outside a partition - either
        % way we apply maxent prior/marginalization
        for j = 1:num_sys_nodes
            %if j is not a denominator but it is an input to i then
            %collapse this dimension
            unidircut = (any(num_nodes_indices(i) == M1) && any(j == M2) && strcmp(bfcut_option,'BRcut')) || ...
                            (any(num_nodes_indices(i) == M2) && any(j == M1) && strcmp(bfcut_option,'FRcut'));
            
            if (~any(j == denom_nodes_indices)||unidircut) && any(j == numerator_nodes(i).input_nodes)
                if any(j == extNodes)
                    past_conditioning_indices = conditioning_indices;
                    past_conditioning_indices{j} = past_state(j) + 1;
                    next_num_node_distribution = next_num_node_distribution(past_conditioning_indices{:});
                else
                    next_num_node_distribution = ...
                       sum(next_num_node_distribution,j)./size(next_num_node_distribution,j);
                end
            end
        end
        
        % the magic
        numerator_conditional_joint = bsxfun(@times,numerator_conditional_joint,next_num_node_distribution);
    end
    
    % conditioning on fixed nodes
    perspective = numerator_conditional_joint ./ sum(numerator_conditional_joint(:));
    
    
    
% P(denom_nodes_f | num_nodes_c = numerator_state) = P(denom_nodes_c | num_nodes_p = numerator_state)
elseif strcmp(bf_option,'forward')
    
    denom_nodes_shift = denom_nodes_indices + num_sys_nodes;
    denom_nodes = nodes(denom_nodes_shift);
    % This is just to define the final size of the distribution
    denom_conditional_joint_size = ones(1,2*num_sys_nodes);
    denom_conditional_joint_size(denom_nodes_indices + num_sys_nodes) = [denom_nodes.num_states];
    denom_conditional_joint = ones(denom_conditional_joint_size);
    denom_inputs = [];
    for i = 1:length(denom_nodes)
        %denom_inputs = union(denom_inputs,denom_nodes(i).input_nodes);
        %denom_inputs = unique([denom_inputs denom_nodes(i).input_nodes]);
        denom_inputs = sort([denom_inputs denom_nodes(i).input_nodes]);
        denom_inputs(denom_inputs((1:end-1)') == denom_inputs((2:end)')) = [];  % faster inline implementation of unique
    end
      
    conditioning_indices = cell(1,2*num_sys_nodes);
    conditioning_indices(:) = {':'};

%     % marginalize over nodes not in numerator, these nodes are outside the
%     % system for this iteration or they are outside a partition - either
%     % way we apply maxent prior/marginalization
    for j = 1:num_sys_nodes
            % condition also on external elements at t0 for op_extNodes == 0 (Freeze), only then is extNodes ~empty.
            if (any(j == num_nodes_indices) || any(j == extNodes)) && (any(j == denom_inputs))
                conditioning_indices{j} = numerator_state(j) + 1;
            end
    end    
    
    for i = 1:length(denom_nodes) %Loop over denominator nodes
        
        next_denom_node_distribution = denom_nodes(i).cpt;
        
        % marginalize over nodes not in denom, these nodes are outside the
        % system for this iteration or they are outside a partition - either
        % way we apply maxent prior/marginalization
        for j = num_sys_nodes:-1:1  %numerator nodes    %Larissa: Why backwards??
            unidircut = (any(denom_nodes_indices(i) == M1) && any(j == M2) && strcmp(bfcut_option,'BRcut')) || ...
                            (any(denom_nodes_indices(i) == M2) && any(j == M1) && strcmp(bfcut_option,'FRcut'));
            % not a numerator or an external node, or it is cut AND it is an input then marginalize             
            if (~(any(j == num_nodes_indices)||any(j == extNodes)) || unidircut) && any(j == denom_nodes(i).input_nodes)
                next_denom_node_distribution = ...
                    sum(next_denom_node_distribution,j)./size(next_denom_node_distribution,j);
                if any(j == num_nodes_indices) && unidircut == 1
                    % average over dimension, but don't collapse dimension as
                    % there is a conditioning index for it (because it's a numerator)
                    next_denom_node_distribution = cat(j, next_denom_node_distribution, next_denom_node_distribution);
                end    
            end
        end
        
        % the magic
        denom_conditional_joint = bsxfun(@times,denom_conditional_joint,next_denom_node_distribution);
    end
    
    
    % conditioning on fixed nodes
    denom_conditional_joint = denom_conditional_joint(conditioning_indices{:});
    permute_order = [num_sys_nodes+1:2*num_sys_nodes 1:num_sys_nodes];
    perspective = permute(denom_conditional_joint,permute_order);

end