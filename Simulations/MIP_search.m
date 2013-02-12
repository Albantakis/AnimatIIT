function [Big_phi_MIP MIP Big_phi_cand MIP_cand] = MIP_search(M,N,Big_phi_M,M_IRR_M,prob_M, phi_M,options)

%%
% Find the Big-phi MIP in a subset M
% M: a subset where Big_phi_MIP is computed
% N: number of elements in the whole system
% Big_phi_M: Big_phi values in every subset, M
% prob - distributions for the concept for each purview
% phi - phi values for each purview in prob

%%

% save('prob_M.mat','prob_M')

% global grain;
% %debug remove me
% fprintf('----------------------------------------------\n');
% fprintf('M = %s\n',mat2str(M));

op_big_phi = options(11);
op_sum = options(12);
op_normalize = options(13);
op_big_phi_dist = options(17);
op_console = options(10);

N_M = length(M);

N_Bp = 0;
for i=1: floor(N_M/2)
    N_Bp = N_Bp + nchoosek(N_M,i);
end

Big_phi_cand = zeros(N_Bp,2);
MIP_cand = cell(N_Bp,1);

whole_i = trans_M(M,N);
Big_phi_w = Big_phi_M(whole_i);

if (any(op_big_phi == [4 5]))
    
    phi_whole = phi_M{whole_i}(:,1)';
    phi_w_concepts = phi_whole(phi_whole ~= 0);
    IRR_whole = M_IRR_M{whole_i};
    concepts_whole_p = zeros(2^N_M,length(phi_w_concepts));
    concepts_whole_f = zeros(2^N_M,length(phi_w_concepts));
    
    z = 1;
    for i = 1:length(phi_whole)
        if (phi_whole(i) ~= 0)
            
            if ~isempty(prob_M{whole_i,1}{i}{1})
                concepts_whole_p(:,z) = prob_M{whole_i,1}{i}{1};
            end
            if ~isempty(prob_M{whole_i,1}{i}{2})
                concepts_whole_f(:,z) = prob_M{whole_i,1}{i}{2};
            end
            z = z + 1;
        end
    end  
    
    whole_concepts_cell_p = cell(length(phi_w_concepts),2);
    whole_concepts_cell_f = cell(length(phi_w_concepts),2);
    
    for i = 1:length(phi_w_concepts)
        
       whole_concepts_cell_p{i,1} = concepts_whole_p(:,i);
       whole_concepts_cell_f{i,1} = concepts_whole_f(:,i); 
       whole_concepts_cell_p{i,2} = phi_w_concepts(i);
       whole_concepts_cell_f{i,2} = phi_w_concepts(i);
       
    end
%     
%     phi_whole = phi_w_concepts;
    
end
    

% are we doing both sides of the partition!?!? <-- YES WHEN IT IS HALF AND
% HALF
l = 1;
for i=1: floor(N_M/2)
    C = nchoosek(M,i);
    N_C = size(C,1);
    for j=1: N_C
        M1 = C(j,:);
        M2 = pick_rest(M,M1);
        
        M1_i = trans_M(M1,N);
        M2_i = trans_M(M2,N);
        
        %debug remove me
%         fprintf('Partition: %s x %s\n',mod_mat2str(M1),mod_mat2str(M2)); 
        
        if(op_sum == 1 || op_big_phi == 0)
            
            Big_phi_partition = Big_phi_M(M1_i) + Big_phi_M(M2_i);
                      
        elseif(op_big_phi == 1)
            
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];

            if (~all(phi == 0))

                concepts = zeros(2^N_M,sum(phi~=0));
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= length(M1_IRR))
                            concepts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        else
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        end
                        z = z + 1;

                    end

                end

%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');
%         disp(concepts);
%         disp(concept_phis);
%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');
%         disp('!!!!!!!!!!!!!');

                Big_phi_partition = big_phi_volume(concepts,concept_phis,grain);

            else
                Big_phi_partition = 0;
            end
            
        elseif (op_big_phi == 2)
            
            M1_IRR = M_IRR_M{M1_i};
            M2_IRR = M_IRR_M{M2_i};

            nIRR = length(M1_IRR) + length(M2_IRR);
            IRRs = cell(nIRR,1);
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            
            for x = 1:nIRR
                if (x <= length(M1_IRR))
                   IRRs{x} = M1_IRR{x};
                else
                   IRRs{x} = M2_IRR{x - length(M1_IRR)};
                end
            end
            

            if (~all(phi == 0))

                concepts = zeros(2^N_M,nIRR);
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= length(M1_IRR))
                            concepts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        else
                            concepts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        end
                        z = z + 1;

                    end

                end

                
                Big_phi_partition = big_phi_info(IRRs,concepts,concept_phis);
                
            else

                Big_phi_partition = 0;
            end
        
        elseif(op_big_phi == 3)
            
            
            
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            nIRR = sum(phi ~= 0);
            
%             disp('***************')
%             disp(M1_IRR);
%             disp(length(M1_IRR));
%             disp(M2_IRR);
%             disp(length(M2_IRR));
%             disp(nIRR)
%             disp(phi)
%             disp(phi ~= 0);
%             disp(sum(phi ~= 0));
            
            if (nIRR > 1)

                concepts_past = zeros(2^N_M,nIRR);
                concepts_future = zeros(2^N_M,nIRR);
                concept_phis = phi(phi ~= 0);

                z = 1;
                for k = 1:length(phi)

                    if (phi(k) ~= 0)
                        
                        if(z <= sum(phi_M{M1_i}(:,1) ~= 0))
                            concepts_past(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                            concepts_future(:,z) = expand_prob(prob_M{M1_i,1}{k}{2},M,M1);
                        else
                            concepts_past(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                            concepts_future(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2},M,M2);
                        end
                        z = z + 1;

                    end

                end
                display = 0;
                Big_phi_partition = big_phi_spacing(concepts_past,concept_phis,display) + big_phi_spacing(concepts_future,concept_phis,display);

            elseif (sum(phi ~= 0) == 1)
                Big_phi_partition = phi((phi ~= 0));
            else
                Big_phi_partition = 0;
            end
            
        elseif (op_big_phi == 4)
            
            M1_IRR = M_IRR_M{M1_i};
            M2_IRR = M_IRR_M{M2_i};

            nIRR = length(M1_IRR) + length(M2_IRR);
            IRR_parts = cell(nIRR,1);
            phi_parts_all = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            
            for x = 1:nIRR
                if (x <= length(M1_IRR))
                   IRR_parts{x} = M1_IRR{x};
                else
                   IRR_parts{x} = M2_IRR{x - length(M1_IRR)};
                end
            end
            

                
            concepts_past = zeros(2^N_M,nIRR);
            concepts_future = zeros(2^N_M,nIRR);
            phi_parts = phi_parts_all(phi_parts_all ~= 0);

            z = 1;
            
            for k = 1:length(phi_parts_all)

                if (phi_parts_all(k) ~= 0)

                    if(z <= sum(phi_M{M1_i}(:,1) ~= 0))
                        if ~isempty(prob_M{M1_i,1}{k}{1})
                            concepts_past(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        end
                        if ~isempty(prob_M{M1_i,1}{k}{2})
                            concepts_future(:,z) = expand_prob(prob_M{M1_i,1}{k}{2},M,M1);
                        end
                    else
                        if ~isempty(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1})
                            concepts_past(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        end
                        if ~isempty(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2})
                            concepts_future(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2},M,M2);
                        end
                    end
                    
                    z = z + 1;

                end

            end

%                 concepts_parts = zeros(2^N_M,nIRR);
%                 phi_parts = phi(phi ~= 0);
% 
%                 z = 1;
%                 for k = 1:length(phi)
% 
%                     if (phi(k) ~= 0)
%                         
%                         if(z <= length(M1_IRR))
%                             concepts_parts(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
%                         else
%                             concepts_parts(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i})}{1},M,M2);
%                         end
%                         z = z + 1;
% 
%                     end
% 
%                 end

%             fprintf('-------------------------------------------------------\n');
%             fprintf('M = %s with partition %s - %s\n',mod_mat2str(M),mod_mat2str(M1),mod_mat2str(M2));
            d_Big_phi = big_phi_shift(M1_IRR, M2_IRR, N, M, IRR_whole,concepts_whole_p,concepts_whole_f,phi_w_concepts, M1, M2,...
                              IRR_parts,concepts_past,concepts_future, prob_M, phi_M{M1_i}(:,1)', phi_M{M2_i}(:,1)', phi_parts,op_big_phi_dist);
                          
        elseif(op_big_phi == 5)
            
            
            
            phi = [phi_M{M1_i}(:,1)' phi_M{M2_i}(:,1)'];
            nIRR = sum(phi ~= 0);
            
%             disp('***************')
%             disp(M1_IRR);
%             disp(length(M1_IRR));
%             disp(M2_IRR);
%             disp(length(M2_IRR));
%             disp(nIRR)
%             disp(phi)
%             disp(phi ~= 0);
%             disp(sum(phi ~= 0));
            


            concepts_past = zeros(2^N_M,nIRR);
            concepts_future = zeros(2^N_M,nIRR);
            concept_phis = phi(phi ~= 0);

            z = 1;
            for k = 1:length(phi)

                if (phi(k) ~= 0)

                    if(z <= sum(phi_M{M1_i}(:,1) ~= 0))
                        concepts_past(:,z) = expand_prob(prob_M{M1_i,1}{k}{1},M,M1);
                        concepts_future(:,z) = expand_prob(prob_M{M1_i,1}{k}{2},M,M1);
                    else
                        concepts_past(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{1},M,M2);
                        concepts_future(:,z) = expand_prob(prob_M{M2_i,1}{k - length(phi_M{M1_i}(:,1))}{2},M,M2);
                    end
                    z = z + 1;

                end

            end
            
            part_concepts_cell_p = cell(nIRR,2);
            part_concepts_cell_f = cell(nIRR,2);

            for k = 1:nIRR

               part_concepts_cell_p{k,1} = concepts_past(:,k);
               part_concepts_cell_f{k,1} = concepts_future(:,k); 
               part_concepts_cell_p{k,2} = concept_phis(k);
               part_concepts_cell_f{k,2} = concept_phis(k);

            end
            
            d_Big_phi = C_distance(part_concepts_cell_p,whole_concepts_cell_p)...
                + C_distance(part_concepts_cell_f,whole_concepts_cell_f);
        
        end
        
        if (any(op_big_phi == [0 1 2 3]))
            
            % 07-06-12 CHANGED TO ABS VALUE
            d_Big_phi = abs(Big_phi_w - Big_phi_partition);
            
        end
        
            
%         
%        Norm = min(length(M1),length(M2)) + min(length(M1),length(M2));
        % Norm = 1; % No normalization
        Norm = 2^min(length(M1),length(M2))-1;
        Big_phi_cand(l,1) = d_Big_phi;
        Big_phi_cand(l,2) = d_Big_phi/Norm;
        MIP_cand{l} = M1;
        
%         if N_M == N
%             fprintf('M1=%s-%s: ',mod_mat2str(M1),mod_mat2str(M2));
%             fprintf('%f-(%f+%f)=%f %f\n',Big_phi_w,Big_phi1,Big_phi2,d_Big_phi,d_Big_phi/Norm);
%         end
        
        l = l + 1;
    end
    
end

% if length(M) == N
%     
%     x = cat(1,concepts_whole_p',concepts_past');
%     nWholeConcepts = size(concepts_whole_p,2);
%     save('sample_partition_4n_sys.mat','x','nWholeConcepts')
% %     figure(1)
% %     conceptscatter(all_concepts,size(concepts_whole_p,2))
%     
% end
    

if (op_normalize == 1 || op_normalize == 2) % Option to normalize or not for new methods of computing big phi
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,2));
else
    [min_norm_Big_phi i_phi_min] = min(Big_phi_cand(:,1));
end

if (op_normalize == 0 || op_normalize == 1)
    Big_phi_MIP = Big_phi_cand(i_phi_min,1);
else
    Big_phi_MIP = Big_phi_cand(i_phi_min,2);
end

MIP = MIP_cand{i_phi_min};
M2 = pick_rest(M,MIP);

if (op_console && Big_phi_MIP ~= 0)
    fprintf('M = %s\nMIP = %s-%s, Big_phi_MIP = %f\n',mat2str(M),mod_mat2str(MIP),mod_mat2str(M2),Big_phi_MIP);
end