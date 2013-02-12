function [H H1 H2] = KLD(prob,prob2)

if any((prob ~= 0) .* (prob2 == 0))
    warning('KL-Divergence may not be correct because prob(i) ~= 0 and prob2(i) == 0 for some i');
end

prob2(prob2==0) = 1; % avoid log0 when computing entropy
H1 = - sum(prob.*log2(prob2)) ;

prob(prob==0) = 1;
H2 = - sum(prob.*log2(prob));
H = H1 - H2;


% prob2 = prob2 + eps; % avoid log0 when computing entropy
% H1 = - sum(prob.*log2(prob2)) ;
% 
% prob = prob + eps;
% H2 = - sum(prob.*log2(prob));
% H = H1 - H2;

