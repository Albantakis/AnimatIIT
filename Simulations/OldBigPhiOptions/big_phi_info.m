function big_phi = big_phi_info(M_IRR,IRR_REP,IRR_phi)

big_phi = 0;

N = size(IRR_REP,1);


% disp(M_IRR)
% disp(IRR_REP)
% disp(IRR_phi)

uniform = ones(N,1)./N;
nIRR = length(M_IRR);

% dist_matrix = gen_dist_matrix(size(IRR_REP,1));
% 
% min_dist = Inf;
% 
% compute pairwise distances
% for i = 1:nIRR
%     for j = i+1:nIRR
%         disp(['FROM' num2str(i) ':'])
%         disp(M_IRR{i});
%         disp(['TO' num2str(j) ':'])
%         disp(M_IRR{j});
%         dist = emd_hat_gd_metric_mex(IRR_REP(:,i),IRR_REP(:,j),dist_matrix);
%         disp(dist);
%         if (dist < min_dist)
%             min_dist = dist;
%             min_i = i;
%             min_j = j;
%         end
%     end
% end
% 
% if (min_dist ~= Inf)
% 
%     disp('SHORTEST DIST FROM:')
%     disp(M_IRR{min_i})
%     disp('TO:')
%     disp(M_IRR{min_j})
%     disp(['= ' num2str(min_dist)])
%     
% end

for i = 1:nIRR
%     disp('THIS PURV');
    this_purview = M_IRR{i};
%     disp('POINT');
    this_concept_point = IRR_REP(:,i);
    
    entropy_loss = 0;
    
    if length(this_purview) == 1
        
        entropy_loss = discrete_entropy(uniform) - discrete_entropy(this_concept_point);
        
    else
        
        entropy_loss_sum = 0;
        count = 0;
        
        for j = 1:nIRR
            
%             ('OTHER PURV')
            other_purview = M_IRR{j};
            
            if ((length(this_purview) - length(other_purview) == 1) &&...
                    all(ismember(other_purview,this_purview)))
                
%                 ('OTHER POINT')
                other_concept_point = IRR_REP(:,j);
                
%                 disp(discrete_entropy(other_concept_point));
%                 disp(discrete_entropy(this_concept_point));
                entropy_loss_sum = entropy_loss_sum + (discrete_entropy(other_concept_point) - discrete_entropy(this_concept_point));
                count = count + 1;
            end
        end
        
        if (count~=0)
            entropy_loss = entropy_loss_sum/count;
        end
   
    end
  
    big_phi = big_phi + entropy_loss * IRR_phi(i);
    
end


end