function prob = partial_prob_comp(x0_b,x1_b,x1,p,b_table,op_fb,M,C_j)
% compute the conditional probability p(x0_b | x1_b)
%
% x1: data, x0_b: partition of x0, x1_b: partition of x1
% p: probability matrix in the whole system



N = log2(size(p,1)); % number of elements in the whole system
% two_pow = zeros(N,1);
% for i=1: N
%     two_pow(i) = 2^(i-1);
% end

two_pow = 2.^(0:(N-1))';

N0_b = length(x0_b); % number of elements in x0
x0_r = 1:N;
x0_r(x0_b) = []; % the rest of x0
N0_r = length(x0_r); % number of elements in the rest of x0

N1_b = length(x1_b);
x1_r = 1:N;
x1_r(x1_b) = []; % the rest of x1
N1_r = length(x1_r); % number of elements in the rest of x1

% x1_i1 = sum(two_pow(x1_b).*x1(x1_b));

if N0_r == 0
    x0_i2_vec = 0;
else
    x0_i2_vec = zeros(2^N0_r,1); % a vector the size of the number of states over the complement of x0_b
    for j=1: 2^N0_r
        % x0_rs = trans2(j-1,N0_r);
        x0_rs = b_table{j,N0_r}; %this gets the column of b_table that holds binary values the size of N0_r
        x0_i2 = sum(two_pow(x0_r).*x0_rs); % this converts back to decimal? - doesn't appear to happen anyway
        x0_i2_vec(j) = x0_i2;
    end
end

if N1_r == 0
    x1_i2_vec = 0;
else
    x1_i2_vec = zeros(2^N1_r,1);
    for j=1: 2^N1_r
        x1_rs = b_table{j,N1_r};
        x1_i2 = sum(two_pow(x1_r).*x1_rs);
        x1_i2_vec(j) = x1_i2;
    end
end

%% original version
if op_fb == 1
    %% backward computation
    % p: 2^N by N
    prob = zeros(2^N0_b,1); % output
    for i = 1:2^N0_b % states of x0 inside the partition
        prob(i) = 1;
        x0_bs = b_table{i,N0_b};
        x0_i1 = sum(two_pow(x0_b).*x0_bs);
        for k = 1:N1_b %% index of x1 neurons
            x1_s = x1(x1_b(k)); % state of the x1 neuron
            p_k = 0;
            for j=1: 2^N0_r % summation over the rest of x0 (outside the partition)
                x0_i2 = x0_i2_vec(j);
                x0_i = 1+x0_i1 + x0_i2;
                if x1_s == 1
                    p_k = p_k + p(x0_i,x1_b(k));
                else
                    p_k = p_k + (1-p(x0_i,x1_b(k)));
                end
            end
            prob(i) = prob(i)*p_k;
        end
    end
else
    %% forward computation
    prob = zeros(2^N1_b,1); % the coditional entropy p(x0_b(not fixed)|x1_b(fixed))
    if op_fb == 0
        % partition
        source = [];
        for j=1: length(C_j)
            source =  [source x0_b(x0_b==C_j(j))];
        end
        x0_in = pick_rest(x0_b,source); % the rest of x0 in M
        x0_out = pick_rest(1:N,x0_b); % the rest of x0 outside M
    elseif op_fb == 2
        % perspective
        source = x0_b;
        x0_in = pick_rest(M,x0_b); % the rest of x0 in M
        x0_out = pick_rest(1:N,M); % the rest of x0 outside M
    end
    N_in = length(x0_in);
    N_out = length(x0_out);
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
    
    % index of the source, Note: interpret x1 as current state  
    if isempty(source) == 1
        x0_so_i = 0;
    elseif length(x1) ~= N
        x0_so_i = sum(two_pow(source).*x1);
    else
        x0_so_i = sum(two_pow(source).*x1(source));
    end
    
    for i=1: 2^N1_b % states of x1 inside the partition
        x1_bs = b_table{i,N1_b};
        prob(i) = 0;
        % fprintf('%s\n',mat2str(x0_bs))
        for j=1: 2^N_in % summation over the rest of x0 in M
            x0_in_i = x0_in_vec(j);
            % fprintf('%d\n',x0_in_i);
            temp = 1;
            for k=1: N1_b % multiplication over x1 states
                x1_s = x1_bs(k);
                p_k = 0;
                for l=1: 2^N_out % summation over the rest of x0 outside M
                    x0_out_i = x0_out_vec(l);
                    x0_i = 1+x0_so_i+x0_in_i + x0_out_i;
                    if x1_s == 1
                        p_k = p_k + p(x0_i,x1_b(k));
                    else
                        p_k = p_k + (1-p(x0_i,x1_b(k)));
                    end
                    % fprintf('k=%d l=%d p_k=%f\n',l,k,p_k);
                end
                temp = temp*p_k;
            end
            % fprintf('temp=%f\n',temp);
            prob(i) = prob(i) + temp;
        end
        % fprintf('prob(i)=%f\n',prob(i));
    end
    
end


%% Normalization
if sum(prob) ~= 0
    prob = prob/sum(prob);
end