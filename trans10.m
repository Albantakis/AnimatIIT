function x = trans10(sigma)

N = length(sigma);

% if nargin < 2    
    two_pow = 2 .^ (0:N-1)';
%     two_pow = two_pow';
%     two_pow = zeros(N,1);
%     for i=1: N
%         two_pow(i) = 2^(i-1);
%     end
% end

% x = 1;
% for i=1: N
%     x = x + 2^(i-1)*sigma(i);
% end

x = 1 + sum(two_pow.*sigma);