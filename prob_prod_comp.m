function prob_prod = prob_prod_comp(prob1,prob2,whole_set,x0_p1,op_fb)

% COMPUTES THE DISTRIBUTION OVER THE STATES OF THE WHOLE SET GIVEN THE DISTRIBUTIONS OF
% THE STATES OF THE SUBSETS
% prob1 = distribution over subset x0_p1
% prob2 = distribution over complement of subset x0_p1
% whole_set = the "whole" system
% x0_p1 = one of the partitions
% op_fb = currently unnecessary 

N = length(whole_set);

if nargin < 5
    op_fb = 0;
end

if isempty(prob2) == 1
    prob_prod = prob1;
elseif isempty(prob1) == 1
    prob_prod = prob2;
else 
    
    N1 = length(x0_p1);
    N2 = N - N1;
    x0_p1_or = x0_p1;
    N_vec = 1:N;
    for i= 1:N1
        x0_p1(i) = N_vec(whole_set==x0_p1_or(i)); % re-index of x0_p1 in terms of whole_set
    end
    x0_p2 = 1:N;
    x0_p2(x0_p1) = [];
    
    if op_fb == 3
        prob1_s = reshape(prob1,[2^N1 2^N1]);
        prob2_s = reshape(prob2,[2^N2 2^N2]);
        prob_prod = zeros(2^N,2^N);
        for i = 1:2^N
            xp_bs = trans2(i-1,N);
            xp_bs1 = xp_bs(x0_p1);
            xp_bs2 = xp_bs(x0_p2);
            xp_i1 = trans10(xp_bs1);
            xp_i2 = trans10(xp_bs2);
            for j=1: 2^N
                xf_bs = trans2(j-1,N);
                xf_bs1 = xf_bs(x0_p1);
                xf_bs2 = xf_bs(x0_p2);
                xf_i1 = trans10(xf_bs1);
                xf_i2 = trans10(xf_bs2);
                prob_prod(i,j) = prob1_s(xp_i1,xf_i1)*prob2_s(xp_i2,xf_i2);
            end
        end
        prob_prod = prob_prod(:);
    else
        prob_prod = zeros(2^N,1);
        x0_bs = zeros(N,1);
        for i=1: 2^N1
            x0_bs1 = trans2(i-1,N1);
            x0_bs(x0_p1) = x0_bs1;
            for j=1: 2^N2
                x0_bs2 = trans2(j-1,N2);
                x0_bs(x0_p2) = x0_bs2;
                x0_i = trans10(x0_bs);
                prob_prod(x0_i) = prob1(i)*prob2(j);
            end
        end
        
    end
end
