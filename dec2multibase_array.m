function multibase_array = dec2multibase_array(decimal_num,state_size_vec)

multibase_array = mod(floor(decimal_num ./ [1 cumprod(state_size_vec(1:end-1))]), state_size_vec);