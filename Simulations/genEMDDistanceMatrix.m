function DistMat = genEMDDistanceMatrix(wholeDists,cutDists, maxEnt) %past whole and cut distributions
%wholeDists and cutDists are matrices of distributions (each column one
%distribution)

wholeDists = [wholeDists; maxEnt];
cutDists = [cutDists; maxEnt];

% should be equal!
rows = size(wholeDists,1);
cols = size(cutDists,1); 

DistMat = zeros(rows, cols); 
for i = 1:rows
    for j = 1:cols
        DistMat(i,j) = emd_hat_gd_metric_mex(wholeDists{i},cutDists{j},gen_dist_matrix(length(wholeDists{i})));

    end    
end
end