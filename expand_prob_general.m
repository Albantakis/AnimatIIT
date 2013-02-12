function joint_distribution = expand_prob_general(prob_in,M,x_in,prob_out)


N = length(M);
N_in = length(x_in); % units inside perspective
N_out = N - N_in;

joint_distribution = zeros(2^N,1); % expanded probability distribution

x_in_re = zeros(N_in,1);
N_vec = 1:N;
for i=1: N_in
    x_in_re(i) = N_vec(M==x_in(i)); % reindex of x in terms of M
end
x_out_re = 1:N;
x_out_re(x_in_re) = [];

x_b = zeros(N,1);

if isempty(prob_in) ~= 1
    for i=1: 2^N_in
        x_b_in = trans2(i-1,N_in);
        x_b(x_in_re) = x_b_in;
        for j=1: 2^N_out
            x_b_out = trans2(j-1,N_out);
            x_b(x_out_re) = x_b_out;
            x0_i = trans10(x_b);
            joint_distribution(x0_i) = prob_in(i)*prob_out(j);
        end
    end
else
    joint_distribution = prob_out;
end

joint_distribution = joint_distribution/sum(joint_distribution);