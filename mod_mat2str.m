function x_s = mod_mat2str(x)
N = length(x);
if N == 0
    x_s = '[]';
elseif N == 1
    x_s = ['[' int2str(x) ']'];
else
    x_s = mat2str(x);
end
