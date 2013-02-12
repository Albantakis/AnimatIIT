function big_phi = big_phi_spacing(concepts,c_phi,display)

nConcepts = size(concepts,2);
nVec = 1:nConcepts;
dim = size(concepts,1);
N = log2(dim);

emd_dist_matrix = gen_dist_matrix(dim);

adj_matrix = zeros(nConcepts);
v_matrix = zeros(nConcepts);

radius = .5;


% compute pairwise distances and if less than .5, subtract appropriate
% amount
for i = 1:nConcepts
    for j = i+1:nConcepts

        dist = emd_hat_gd_metric_mex(concepts(:,i),concepts(:,j),emd_dist_matrix);
       
        if (dist < radius)
            
            v = 2 * (radius - dist);
            adj_matrix(i,j) = 1;
            adj_matrix(j,i) = 1;
            v_matrix(i,j) = v;
            v_matrix(j,i) = v;
            
        end
    end
end

% beginning of code to take out redundant concepts
% [I,J] = find(v_matrix == 1);
% for i = 1:length(I)
%     [min,index] = 

cliques = logical(maximalCliques(adj_matrix));

if (display)
disp(v_matrix);
disp(cliques);
disp(concepts');
disp(c_phi);end

max_clique_sizes = sum(cliques,1);
max_clique_size_total = max(max_clique_sizes);

clique_sets = cell(max_clique_size_total,1);

for i = 1:size(cliques,2)
    
    for j = 2:max_clique_sizes(i)
        
        concepts = nVec(cliques(:,i));
        clique_sets{j} = union_edit(clique_sets{j}, nchoosek(concepts,j), 'rows');
        
    end
end

% big_pi = zeros(max_clique_size_total,1);
big_phi = sum(c_phi);
% big_phi = nConcepts;
%for cliques of each size
for k = 2:max_clique_size_total
    
    big_pi = 0;
    
    these_cliques = clique_sets{k};
    
    %for each clique of this size
    for m = 1:size(these_cliques,1)
        
        this_clique = these_cliques(m,:);
        product = 1;
        for i = 1:k
            for j = i+1:k
                
                row = this_clique(i);
                col = this_clique(j);
                product = product * v_matrix(row,col);
                
            end
        end
        
        if (display) disp(product); end
        
%         big_pi = big_pi + product^(1 / nchoosek(k,2));
        big_pi = big_pi + product^(1 / nchoosek(k,2)) * min(c_phi(this_clique));

    end
    if (display)
    disp(k);
    disp(product);
    disp(big_pi);
    end
    big_phi = big_phi + (-1)^(k-1) * big_pi;
end

% big_phi = sum(c_phi) + sum(big_pi .* (ones(max_clique_size_total,1) * (-1)).^(1:max_clique_size_total));

        
        









% THE COMMMENTED OUT CODE IS AN OLD WAY
% overlap = 0;
% 
% % compute pairwise distances and if less than .5, subtract appropriate
% % amount
% for i = 1:nConcepts
%     for j = i+1:nConcepts
% 
%         dist = emd_hat_gd_metric_mex(concepts(:,i),concepts(:,j),emd_dist_matrix);
%        
%         if (dist < .5)
%             
%             overlap = overlap + (.5 - dist) * min(c_phi(i),c_phi(j));
%             
%         end
%     end
% end
% 
% scaling_factor = nConcepts / nchoosek(nConcepts,2);
% 
% big_phi = nConcepts - overlap * scaling_factor;
% 
% end

