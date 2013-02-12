function simplex_points = simplex_volume(nElements,steps_per_dimension)

% C = N+1 x Nc matrix of concepts as columns, N = dimension of simplex, Nc
%   is the number of concepts

% C_phi = 1xNc vector of phi values for each concept

% steps_per_dimension is how many points are seached in each dimension

% big_phi = BIG PHI!!
% matlabpool close force local
% matlabpool

tic

% get dimension of simplex
N = 2^nElements - 1;


%% Build R

% initialize R, a matrix of coordinates of the simplex in the projected
% space where each column is the coordinates of one of the vertices of the
% simplex

R = zeros(N,N+1);
R(1,1) = 1;
R(1,2:end) = -1/N;

for j = 2:N
    
    R(j,j) = nthroot(sum(R(:,j).^N),N);
    R(j,j+1:end) = -1/N - R(1:j-1,j)'*R(1:j-1,j+1);

end

R(N,N+1) = -R(N,N);

%% Build the transform matrix T which maps from concept space to the space of R

r0_mat = repmat(R(:,1),1,N);
T = inv(R(:,2:end) - r0_mat);



%% Test all points in lattice bounded by [-1,1] for all dimensions with granularity passed in arguments

% numPoints = (2/grain + 1)^N;
% test_points = cell(numPoints,1);

ref_point = ones(N,1) * -1;
distance = zeros(N,1);

nSimplex = 0;
step_size = 2 / (steps_per_dimension - 1);

simplex_points = cell(steps_per_dimension^N,1);

for i = 0:(steps_per_dimension^N-1)
    
    eval(['simplex_points_' num2str(i) ' = cell(steps_per_dimension^N,1);']);
    eval(['count_' num2str(i) ' = 0;']);
    
end

parfor i = 0:(steps_per_dimension^N-1)
    
   for j = 1:N
       
       distance(j) = mod(i,steps_per_dimension);
       i = (i - distance(j)) / steps_per_dimension;
      
   end
       
   test_point = ref_point + distance * step_size;

   r0 = zeros(N,1); r0(1) = 1;

   gammaTemp = T * (test_point - r0);
   gamma = zeros(N+1,1);
   gamma(1) = 1 - sum(gammaTemp);
   gamma(2:end) = gammaTemp;
   
   
   if (sum(gamma) == 1 && all(gamma >= 0))
       
       nSimplex = nSimplex + 1;
%        eval(['simplex_points_' num2str(i) '{count_' num2str(i) '} = gamma;']);
%        eval(['count_' num2str(i) ' = count_' num2str(i) ' + 1;'])
       
       
   end
       
end


toc
% matlabpool close
end