function [big_phi_mip distance_sum phi_sum] = big_phi_shift(M1_IRR, M2_IRR, N, M, IRR_whole,concepts_whole_p,concepts_whole_f,phi_whole,...
    part1, part2, IRR_parts,concepts_parts_p,concepts_parts_f, all_distributions, phi_part1, phi_part2, concept_phi_parts,op_big_phi_dist)

nWholeConcepts = length(IRR_whole);
% nPartConcepts = length(concept_phi_parts);
% partitionedCheck = zeros(nPartConcepts,1);

dist_matrix = gen_dist_matrix(size(concepts_whole_p,1));

distance_sum = 0;
phi_sum = 0;

part1_index = trans_M(part1,N);
part2_index = trans_M(part2,N);

N1 = length(part1);
subsets_part1 = cell(2^N1-1,1);
k = 1;
for i = 1:N1 
    C = nchoosek(part1,i); % create a matrix of combinations of M of size i
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j = 1:N_C % for all combos of size i
        x0 = C(j,:); % pick a combination
        subsets_part1{k} = x0;% store combo
        k = k + 1;
    end
end

N2 = length(part2);
subsets_part2 = cell(2^N2-1,1);
k = 1;
for i = 1:N2 
    C = nchoosek(part2,i); % create a matrix of combinations of M of size i
    N_C = size(C,1);
    % fprintf('i=%d N_c=%d\n',i,N_C);
    for j = 1:N_C % for all combos of size i
        x0 = C(j,:); % pick a combination
        subsets_part2{k} = x0;% store combo
        k = k + 1;
    end
end



for i = 1:nWholeConcepts
    
    
    IRR_w = IRR_whole{i};
    
%     fprintf('Checking Concept: %s\n',mod_mat2str(IRR_w));
    
    % contained in part1
    if all(ismember(IRR_w,part1))
        
        for j = 1:length(subsets_part1)
            
            part1_set = subsets_part1{j};
            
            if (length(IRR_w) == length(part1_set) && all(ismember(IRR_w,part1_set)))
                
                concept_index = j;
                break
            end
        end
        
        
        concept_past = expand_prob(all_distributions{part1_index,1}{concept_index}{1},M,part1);
        concept_future = expand_prob(all_distributions{part1_index,1}{concept_index}{2},M,part1);
        part_phi = phi_part1(concept_index);
        
    % contained in part2    
    elseif all(ismember(IRR_w,part2))
        
        for j = 1:length(subsets_part2)
            
            part2_set = subsets_part2{j};
            
            if (length(IRR_w) == length(part2_set) && all(ismember(IRR_w,part2_set)))
                
                concept_index = j;
                break
            end
        end
%         
%         disp(M)
%         disp(part2)
%         disp(part2_set)
        concept_past = expand_prob(all_distributions{part2_index,1}{concept_index}{1},M,part2);
        concept_future = expand_prob(all_distributions{part2_index,1}{concept_index}{2},M,part2);
        part_phi = phi_part2(concept_index);
    
	% this purview was split by partition
    else
        
        part1_count = sum(ismember(IRR_w,part1));
        for j = 1:length(subsets_part1)
            
            part1_set = subsets_part1{j};
            
            if (part1_count == length(part1_set) && all(ismember(part1_set,IRR_w)))
                
                part1_concept_index = j;
                break
            end
        end
        
        part2_count = sum(ismember(IRR_w,part2));
        for j = 1:length(subsets_part2)
            
            part2_set = subsets_part2{j};
            
            if (part2_count == length(part2_set) && all(ismember(part2_set,IRR_w)))
                
                part2_concept_index = j;
                break
            end
        end
        
        concepts_past_p1 = all_distributions{part1_index,1}{part1_concept_index}{1};
        concepts_future_p1 = all_distributions{part1_index,1}{part1_concept_index}{2};
        
        concepts_past_p2 = all_distributions{part2_index,1}{part2_concept_index}{1};
        concepts_future_p2 = all_distributions{part2_index,1}{part2_concept_index}{2};
        
        concept_past = prob_prod_comp(concepts_past_p1, concepts_past_p2, M, part1, 0);
        concept_future = prob_prod_comp(concepts_future_p1, concepts_future_p2, M, part1, 0);
        
        part_phi = 0;
        
    end
       
        
    if op_big_phi_dist == 0
        
        %add in distances b/w past concepts for this purview
        past_dist = KLD(concepts_whole_p(:,i), concept_past);
        distance_sum = distance_sum + past_dist;
        %add in distances b/w future concepts for this purview
        future_dist = KLD(concepts_whole_f(:,i), concept_future);
        distance_sum = distance_sum + future_dist;

    elseif op_big_phi_dist == 1
        
%         disp(concepts_whole_p(:,i))
%         disp(concept_past)
        %add in distances b/w past concepts for this purview
        past_dist = emd_hat_gd_metric_mex(concepts_whole_p(:,i), concept_past, dist_matrix);
        distance_sum = distance_sum + past_dist;
        %add in distances b/w future concepts for this purview
        future_dist = emd_hat_gd_metric_mex(concepts_whole_f(:,i), concept_future, dist_matrix);
        distance_sum = distance_sum + future_dist;
    end
    
    phi_sum = phi_sum + abs(phi_whole(i) - part_phi);
    
    %for deubbing, take out
%     fprintf('\tDistance to partitioned past distribution: %f\n',past_dist);
%     fprintf('\tDistance to partitioned future distribution: %f\n',future_dist);
%     fprintf('\tWhole Small Phi: %f | Partitioned Small Phi: %f\n',phi_whole(i),part_phi);
%     fprintf('\tSmall Phi Diff: %f(abs) %f(whole - part)\n',abs(phi_whole(i) - part_phi),phi_whole(i) - part_phi);
    
end
        
    % For debuggin take out
%     fprintf('Concept From Whole: %s\n',mod_mat2str(IRR_w));
%     concept_w = concepts_whole(:,i);
%     phi_w = phi_whole(i);
    
%     partner_found = 0;
%     
%     for j = 1:nPartConcepts
%         
%         % we haven't already found this guy's partner
%         if (partitionedCheck(j) == 0)
%             
%             IRR_p = IRR_parts{j};
%         
%             if (length(IRR_p) == length(IRR_w) && all(ismember(IRR_w,IRR_p)))
%                 
%                 partner_found = 1;
%                 partitionedCheck(j) = 1;
%                 
%                 if op_big_phi_dist == 0
%                     %add in distances b/w past concepts for this purview
%                     past_dist = KLD(concepts_whole_p(:,i), concepts_parts_p(:,j));
%                     distance_sum = distance_sum + past_dist;
%                     %add in distances b/w future concepts for this purview
%                     future_dist = KLD(concepts_whole_f(:,i), concepts_parts_f(:,j));
%                     distance_sum = distance_sum + future_dist;
%                     
%                 elseif op_big_phi_dist == 1
%                     %add in distances b/w past concepts for this purview
%                     past_dist = emd_hat_gd_metric_mex(concepts_whole_p(:,i), concepts_parts_p(:,j),dist_matrix);
%                     distance_sum = distance_sum + past_dist;
%                     %add in distances b/w future concepts for this purview
%                     future_dist = emd_hat_gd_metric_mex(concepts_whole_f(:,i), concepts_parts_f(:,j),dist_matrix);
%                     distance_sum = distance_sum + future_dist;
%                 
%                 end
%                 
% %                 %for deubbing, take out
% %                 fprintf('\tDistance to past distribution: %f\n',past_dist);
% %                 fprintf('\tDistance to future distribution: %f\n',future_dist);
% %                 fprintf('\tSmall Phi Diff: %f(abs) %f(whole - part)\n',abs(phi_whole(i) - concept_phi_parts(j)),phi_whole(i) - concept_phi_parts(j));
%                 %add in phi difference
%                 phi_sum = phi_sum + abs(phi_whole(i) - concept_phi_parts(j));
%                 
%             end
%         end
%         
%     end
%     
%     % if we didn't find a partner, just add in the phi value
%     if ~partner_found
% %         fprintf('\tConcept does not exist in partitioned system\n');
% %         fprintf('\tSmall Phi Contribution: %f\n',phi_whole(i));
%         phi_sum = phi_sum + phi_whole(i);
% 
%     end
% 
%     % if this concept is from 
% end

% for concepts in the partitioned system which do not exist in the whole,
% add their small_phi in
% if (any(partitionedCheck == 0))
%     fprintf('\tPartitioned concepts not in whole:\n');
% end

% for i = 1:nPartConcepts
%     
%     if (partitionedCheck(i) == 0)
%         
%         fprintf('\t\t%s: small phi contribution, %f',mod_mat2str(IRR_parts{i}),concept_phi_parts(i));
%         phi_sum = phi_sum + concept_phi_parts(i);
%         
%     end
% end

big_phi_mip = distance_sum + phi_sum;
                
                

        

