function prob = partial_prob_forward(x0_so,x0_in,x0_out,x1_b,x0,x1,p,b_table)

% compute the conditional probability of x1 given x0, p(x1(fixed)|x0(fixed))
% x0_so: fixed state, x0_in: injecting noise in units (maxEnt)
% x0_out: injecting noise in connections (complete noise)
% x0, x1: given data
% p: transition probability matrix (TPM)
% b_table: table used for converting binary sequences into decimal number

N = size(p,2); % number of elements in the whole system
% two_pow = zeros(N,1);
% for i=1: N
%     two_pow(i) = 2^(i-1);
% end
two_pow = 2.^(0:N-1)';

% this is always 0????
N_in = length(x0_in); % number of elements of x0 inside the partition (maxEnt)
N_out = length(x0_out); % number of elements of x0 outside the partition (complete noise)
N1_b = length(x1_b); % number of elements of x1 in target

x0_in_vec = zeros(2^N_in,1);
x0_out_vec = zeros(2^N_out,1);
if N_in ~= 0
    for j=1: 2^N_in
        x0_s_in = b_table{j,N_in};
        x0_in_i = sum(two_pow(x0_in).*x0_s_in);
        x0_in_vec(j) = x0_in_i;
    end
else
    x0_in_vec = 0;
end

if N_out ~= 0
    for j=1: 2^N_out
        x0_s_out = b_table{j,N_out};
        x0_out_i = sum(two_pow(x0_out).*x0_s_out);
        x0_out_vec(j) = x0_out_i;
    end
else
    x0_out_vec = 0;
end

% index of the denominator
if isempty(x0_so) == 1
    x0_so_i = 0;
else
    x0_so_i = sum(two_pow(x0_so).*x0);
end

prob = 0; % the coditional entropy p(x1(fixed)|x0(fixed))

% this looks like OR(AND(OR)) which is to say SUM OF PRODUCTS OF SUMS
% for each state of x0_in which is always empty set?
for j=1: 2^N_in % summation over the rest of x0 inside the partition
    x0_in_i = x0_in_vec(j); % always 0?
    temp = 1;
    
    % for each element in numerator
    for k = 1:N1_b
        x1_s = x1(k); % state of this element
        p_k = 0; % set OR prob to 0
        for l = 1:2^N_out % summation over the rest of denominator outside the partition
            x0_out_i = x0_out_vec(l); % get the state of rest of denom
            % fprintf('%d %d %d\n',x0_so_i,x0_in_i,x0_out_i);
            % x0_i
            x0_i = 1+x0_so_i+x0_in_i + x0_out_i; % get full state number
            if x1_s == 1
                % add in prob
                p_k = p_k + p(x0_i,x1_b(k));
            else
                p_k = p_k + (1-p(x0_i,x1_b(k)));
            end
        end
        temp = temp*p_k;
    end
    prob = prob + temp;
    % fprintf('j=%d prob=%f temp=%f\n',j,prob,temp)
end

% prob = prob/2^(N_in+N_out);