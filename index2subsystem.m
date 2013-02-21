function subsystem = index2subsystem(index)

number_vector = 1:ceil(log2(index));

bin_string = dec2bin(index,length(number_vector));
bin_logical = false(length(number_vector),1);

for i = 1:length(number_vector)
    
    bin_logical(i) = logical(str2num(bin_string(i)));
    
end

bin_logical = flipud(bin_logical);

subsystem = number_vector(bin_logical);

if length(subsystem) == 1
    subsystem = subsystem + 1;
end

if isempty(subsystem)
    subsystem = 1;
end