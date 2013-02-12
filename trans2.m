function [binary] = trans2(value,length)
% TRANS2(VALUE,LENGTH) returns an array of LENGTH 0's and 1's, representing 
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
binary = zeros(length,1);
for i = 1:length
    binary(i) = value - 2 * floor(value/2);
    value = floor(value/2);
end
% binary
% toc
