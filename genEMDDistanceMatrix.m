function [Dist, indD] = genEMDDistanceMatrix(Dists, maxEnt, dist_mat) %past whole and cut distributions
%wholeDists and cutDists are matrices of distributions (each column one
%distribution)

for pf = 1:2
    wholeDists = [Dists(:,pf,1); maxEnt(:,pf)];
    cutDists = [Dists(:,pf,2); maxEnt(:,pf)];

    % should be equal!
    rows = size(wholeDists,1);
    cols = size(cutDists,1); 
    if pf == 1
        DistMat = zeros(rows, cols,2); 
    end    
    for i = 1:rows
        for j = 1:cols
            DistMat(i,j,pf) = emd_hat_gd_metric_mex(wholeDists{i},cutDists{j},dist_mat(1:length(wholeDists{i}),1:length(wholeDists{i})));
        end    
    end
end
Dist =  sum(DistMat,3); 
%Larissa: The following does not give much increase in computation speed afterall 
indD = diag(Dist,0) ~= 0;
indD(end) = 1;
Dist = Dist(indD,indD);
end