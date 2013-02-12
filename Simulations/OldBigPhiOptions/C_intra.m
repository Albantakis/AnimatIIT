function [intra] = C_intra(C)

%% compute Intradistance for a constellation
% C: Constellation (cell with first column the distributions, secondo
% column the small phi value)

% option:
if isempty(C)
    intra = 0;
    return
end

opt_EMD = 0;    %1: Use EMD to compute distance between distribution,
                %0: Use the interdistance between distributions
opt_fi_count = 0;   %1: the difference in phi counts for the intra-distance

if opt_fi_count == 1
    Nc=size(C,1);
    N=length(C{1,1});
    tot = 0;
%    Wc=0;
    if opt_EMD==1
        for i=1:Nc
%           Wc=Wc+C{i,2};
           for j =(i+1) : Nc
                   tot = tot + (emd_hat_gd_metric_mex(C{i,1},C{j,1},gen_dist_matrix(N))+abs(C{i,2}-C{j,2})) * sqrt(C{i,2}*C{j,2});
           end
        end
 %       tot = tot/Wc^2;
    else
        for i=1:Nc
           Wc=Wc+C{i,2};
           for j =(i+1) : Nc
                   tot = tot +  k_inter(C{i,1},C{j,1})* C{i,2}*C{j,2};
           end
        end
 %       tot = tot/Wc^2;
    end
    intra = tot;
else
    Nc=size(C,1);
    N=length(C{1,1});
    tot = 0;

    if opt_EMD==1
        for i=1:Nc
           for j =(i+1) : Nc
                   tot = tot +  emd_hat_gd_metric_mex(C{i,1},C{j,1},gen_dist_matrix(N)) * sqrt(C{i,2}*C{j,2});
           end
        end
        if Nc>1 
%            tot = tot/sqrt(nchoosek(Nc,2));
        end
    else
        for i=1:Nc
           for j =(i+1) : Nc
                   tot = tot +  k_inter(C{i,1},C{j,1})* sqrt(C{i,2}*C{j,2});
           end
        end
        if Nc>1 
 %           tot = tot/nchoosek(Nc,2);
        end
    end
    intra = tot;
end
end