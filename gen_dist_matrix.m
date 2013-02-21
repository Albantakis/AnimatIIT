function D = gen_dist_matrix(N)

% generates the distance matrix between vertices which can be used as
% weights we using the Earth Mover's Distance between distributions

% N = support size of simplex, aka it's the number of states in our system,
% aka it's 2^K for a K node system

D = zeros(N,N);


% TODO: I bet I can do this vectorized...
for i = 0:N-1
    for j = 0:N-1
        D(i+1,j+1) = sum(dec2bin(bitxor(i,j)) == '1');
    end
end

