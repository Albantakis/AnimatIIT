function [inter] = C_inter(C,D)

%% compute Interdistace between two constellations using the kernel. 
% This is to be used 3 times to compute the actual "distance" (see function
% C_distance)
%
% P,Q: Constellations (cells with distributions, small_phis)

if isempty(C) || isempty(D)
    inter = 0;
    return
end

opt_EMD = 0;    %1: Use EMD to compute distance between distribution,
                %0: Use the interdistance between distributions

opt_weighted = 1; % 0: not weighted 1: we rescale the sum of the weights to 1 2: Rescaled and take square root  

%%
Nc=size(C,1);
Nd=size(D,1);

N = length(C{1,1});  %Assuming they are all the same length
tot = 0;

if opt_weighted == 0 && opt_EMD==1 
    for i=1:Nc
        for j =1:Nd
            tot = tot + 2^-(emd_hat_gd_metric_mex(C{i,1},D{j,1},gen_dist_matrix(N))+abs(C{i,2}-D{j,2}));
        end
    end
elseif  opt_weighted == 0 && opt_EMD==0 
    for i=1:Nc
        for j =1:Nd
            tot = tot + 2^-(k_distance(C{i,1},D{j,1})+abs(C{i,2}-D{j,2}));
        end
    end
elseif  opt_weighted ~= 0
    Wc = 0;
    for i=1:Nc 
        Wc = Wc + C{i,2}; 
    end
    Wd = 0;
    for i=1:Nd 
        Wd = Wd + D{i,2}; 
    end
    if opt_weighted == 1 && opt_EMD==1 
        for i=1:Nc
            for j =1:Nd
                tot = tot + 2^-(emd_hat_gd_metric_mex(C{i,1},D{j,1},gen_dist_matrix(N))+abs(C{i,2}-D{j,2})) * C{i,2}*D{j,2} ;
            end
        end
        tot = tot/(Wc*Wd);
    elseif  opt_weighted == 1 && opt_EMD==0 
        for i=1:Nc
            for j =1:Nd
                tot = tot + 2^-(k_distance(C{i,1},D{j,1})+abs(C{i,2}-D{j,2})) * C{i,2}*D{j,2};
            end
        end
        tot = tot/(Wc*Wd);        
    elseif opt_weighted == 2 && opt_EMD==1 
        for i=1:Nc
            for j =1:Nd
                tot = tot + 2^-(emd_hat_gd_metric_mex(C{i,1},D{j,1},gen_dist_matrix(N))+abs(C{i,2}-D{j,2})) * sqrt(C{i,2}*D{j,2}) ;
            end
        end
        tot = tot/sqrt((Wc*Wd));               
    elseif  opt_weighted == 2 && opt_EMD==0 
        for i=1:Nc
            for j =1:Nd
                tot = tot + 2^-(k_distance(C{i,1},D{j,1})+abs(C{i,2}-D{j,2})) * sqrt(C{i,2}*D{j,2});
            end
        end
        tot = tot/sqrt(Wc*Wd);               
    end
end
    inter = tot;
end
