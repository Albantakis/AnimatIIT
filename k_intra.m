function [intra] = k_intra(P)

%% compute Intradistance for a distribution
% P: Distribution

%%
op_square_root = 0;  %0: no square root; 1: square root, divide by N

N=length(P);

% if N <= 5
    
%         tic
    A = repmat((0:N-1)',1,N);
    B = A';
    hamming_distance_factors = sum(dec2bin(bitxor(A,B)) == '1',2);

    distribution_products = P * P';
    distribution_products = distribution_products(:);

    intra = sum(distribution_products .* hamming_distance_factors)/2;

%     toc

% else
% tic
%     tot = 0;
% 
%     if op_square_root==1
%         for i=1:N
%            for j =(i+1) : N
%                    tot = tot + sum(dec2bin(bitxor(i-1,j-1)) == '1') * sqrt(P(i)*P(j));
%            end
%         end
%         intra = tot/N;
%     else
%        for i=1:N
%            for j =(i+1) : N
%                     tot = tot + sum(dec2bin(bitxor(i-1,j-1)) == '1') * P(i)*P(j);
%            end
%        end
%        intra = tot
%     end
%     toc
end
% end