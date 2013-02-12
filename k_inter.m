function [inter] = k_inter(P,Q)

%% compute Interdistace between two distributions using the kernel. 
% This is to be used 3 times to compute the actual "distance" (see function
% k_distance)
%
% P,Q: Distributions as column vectors

op_sim = 0;    %0: similarity is 2^-(distance), 1: similarity is 1/(distance+1)
op_sqr = 0;    %0: no square root of prod. of dist. 1: Sqrt of prod.

%% 

N=length(P);




% if (N <= 5)
    
%     tic
    A = repmat((0:N-1)',1,N);
    B = A';
    hamming_distance_factors = 2.^(-sum(dec2bin(bitxor(A,B)) == '1',2));

    distribution_products = P * Q';
    distribution_products = distribution_products(:);

    inter = sum(distribution_products.*hamming_distance_factors);

%     toc
% 
%     tic
%     A = repmat((0:N-1)',1,N);
%     B = A';
%     hamming_distance_factors = exp((-sum(dec2bin(bitxor(A,B)) == '1',2))*log(2));
% 
%     distribution_products = P * Q';
%     distribution_products = distribution_products(:);
% 
%     inter = sum(distribution_products.*hamming_distance_factors)
% 
%     toc
% 
% % else
%     tic
% 
%     tot = 0;
% 
%     if (op_sim==0) && (op_sqr==0)
%         for i=1:N
%             for j =1:N
%                 tot = tot + 2^-(sum(dec2bin(bitxor(i-1,j-1)) == '1')) * P(i)*Q(j);
%             end
%         end
%     elseif (op_sim==0) && (op_sqr==1)
%         for i=1:N
%             for j =1:N
%                 tot = tot + 2^-(sum(dec2bin(bitxor(i-1,j-1)) == '1')) * sqrt(P(i)*Q(j));
%             end
%         end
%     elseif (op_sim==1) && (op_sqr==0)
%         for i=1:N
%             for j =1:N
%                 tot = tot + 1/(1+(sum(dec2bin(bitxor(i-1,j-1)) == '1'))) * P(i)*Q(j);
%             end
%         end
%     elseif (op_sim==1) && (op_sqr==1)
%         for i=1:N
%             for j =1:N
%                 tot = tot + 1/(1+(sum(dec2bin(bitxor(i-1,j-1)) == '1'))) * sqrt(P(i)*Q(j));
%             end
%         end
%     end
%     inter = tot
%     toc
% end
    
end

