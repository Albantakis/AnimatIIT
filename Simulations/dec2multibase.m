function [multidec_array] = dec2multibase(value,state_size_vec)
% TRANS2(VALUE,STATE_SIZE_VEC) returns an array of LENGTH 0's and 1's, representing 
%   the binary number equal to VALUE
%
%   This function also comes as compiled MEX files, the m-file here is a
%   backup if the current 

% SHOULD WE FLIPLR RESULT

% this is a case where the for loop seems to be faster, see for yourself...
% tic
% binary = mod(floor(value./(2.^(0:length-1))),2)
% toc
% 
% tic
num_nodes = length(state_size_vec);
multidec_array = zeros(num_nodes,1);
for i = 1:num_nodes
%     multidec_array(i) = value - 2 * floor(value/2);
    multidec_array(i) = mod(value,state_size_vec(i));
    value = floor(value/state_size_vec(i));
end
% binary
% toc
