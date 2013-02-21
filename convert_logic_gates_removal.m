function new_logic_gates = convert_logic_gates_removal(prev_logic_gates, old_inputs, new_inputs)
% old_inputs and new_inputs are cell arrays
new_logic_gates = prev_logic_gates;


for i = 1:numel(prev_logic_gates)
    if numel(new_inputs{i}) < numel(old_inputs{i}) %if not, nothing changes
        if isempty(new_inputs{i})
            new_logic_gates(i) = 0; % doesn't exist any more, outputs 0 always
        else
            
            logic_type = prev_logic_gates(i);

            % AND
            if logic_type == 1
                    new_logic_gates(i) = 0; % doesn't exist any more
            % OR  elseif logic_type == 2    stays OR, e.g. output = any(input); for one input this is equivalent to COPY
                
            % XOR elseif logic_type == 3    stays XOR, e.g. output = sum(input) == 1; for one input this is equivalent to COPY

            % COPY if no input, it's lost anyways   
            
            % NOT - should have only one input --> is either destroyed or stays the same
            
            % NULL - doesn't have inputs anyways
            
            % MAJORITY -> converted to linear threshold unit with Threshold as before
            elseif logic_type == 7  

                Thres = ceil(numel(old_inputs{i})/2);
                new_logic_gates(i) = 10+Thres;

            % MINORITY -> converted to below threshold unit with Threshold as before
            elseif logic_type == 8
                
                Thres = ceil(numel(old_inputs{i})/2);  %it's still ceil, because the below threshold unit will be "<" not "<="
                new_logic_gates(i) = 20+Thres;

            %PARITY -> gets destroyed (easiest)
            elseif logic_type == 9
                new_logic_gates(i) = 0; % doesn't exist any more
            
            %Linear Threshold units keep their Threshold, i.e. don't change

            end
        end
    end  
end