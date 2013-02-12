function prob = patial_prob_comp_bf(T,x0_so,x0_in,x0_out,x1_in,x1_out,source,p,b_table)

% computing the conditional probability of the past xp and the future xf
% given current state x0, p(xp,xf|x0) 
% x1_b: partition in the past and the future
% x0_so: target in the current state
% x0_in: x0 inside the partition in the current state
% x0_out: x0 outside the partition in the current state
% source: given data
% p: transition probability matrix (TPM)

N1 = length(x1_in);

x0_all = [x0_so x0_in];
x0_all = sort(x0_all);
N0_in = length(x0_in);

x0_state = source;

%% backward computation p(xp|x0)
p_b = zeros(2^N1,2^N0_in);
p_f = zeros(2^N1,2^N0_in);

for i=1: 2^N0_in
    if N0_in ~= 0
        x0_state(x0_in) = b_table{i,N0_in};
    end
    % backward without normalization
    op_fb = 3;
    xp_so = x1_in;
    xp_in = [];
    xp_out = x1_out;
    xc = x0_all; % current state
    p_b(:,i) = partial_prob_comp_time(T,xp_so,xp_in,xp_out,xc,x0_state,p,b_table,op_fb);
    
    % forward with normalization
    op_fb = 0;
    xc_so = x0_all;
    xc_in = [];
    xc_out = x0_out;
    xf = x1_in; % future state
    p_f(:,i) = partial_prob_comp_time(T,xc_so,xc_in,xc_out,xf,x0_state,p,b_table,op_fb);
    
    % fprintf('x0=%s p_b=%s p_f=%s\n',mat2str(x0_state),mat2str(p_b(:,i)),mat2str(p_f(:,i)));
end


prob = p_b*p_f';
Norm = sum(sum(prob));


if Norm ~= 0
    prob = prob/Norm;
end
