function index = state2index(state_vec, num_states_vec)



index = state_vec(1) + 1;

for i = 2:length(num_states_vec)
    
    
   index = index + state_vec(i)*prod(num_states_vec(1:i-1));
    
    
end

    

