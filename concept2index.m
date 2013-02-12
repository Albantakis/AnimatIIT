function index = concept2index(numerator, subsystem)

% size of subsystem
N = length(subsystem);
% size of numerator
M = length(numerator);

% init
index = 0;
% count concepts with smaller numerators
for i = 1:M-1
    
    index = index + nchoosek(N,i);
    
end

% build set of concepts with same size numerator
concepts_size_M = nchoosek(subsystem,M);
% find the concept of interest
[location] = ismember(concepts_size_M,numerator,'rows');

% add to total
index = index + find(location == 1);