function numerator = index2concept(index, subsystem)

% INCOMPLETE FUNCTION 

N = length(subsystem);

M_cell = cell(2^N - 1,1);
    
k = 1;
for i = 1:N 
    C = nchoosek(subsystem,i); 
    N_C = size(C,1);
    for j = 1:N_C % for all combos of size i
        x0 = C(j,:); % pick a combination
        M_cell{k} = x0;% store combo
        k = k + 1;
    end
end

numerator = M_cell{index};
end