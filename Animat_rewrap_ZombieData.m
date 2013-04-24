function data = Animat_rewrap_ZombieData(full_sys,Big_phi_M,phi_M, concept_MIP_M, purviews_M)

data.Phi = Big_phi_M;
num_concepts = 2^length(full_sys) - 1;
% initialize size of concept struct array
data.concept(num_concepts).phi.min = 0;
for concept_i = 1:num_concepts
% for now we set concept to reducible and correct after for-loop
% (see below)
data.concept(concept_i).is_irreducible = 0;
if ~isempty(phi_M)
    % phi values for each concept
    data.concept(concept_i).phi.min = phi_M(concept_i,1);
    data.concept(concept_i).phi.backwards = phi_M(concept_i,2);
    data.concept(concept_i).phi.forwards = phi_M(concept_i,3);
    % denominator
    data.concept(concept_i).denominator.backwards = sort([concept_MIP_M{concept_i}{:,1,1}]);
    data.concept(concept_i).denominator.forwards = sort([concept_MIP_M{concept_i}{:,1,2}]);
end
end
% label irreducible concepts
if ~isempty(purviews_M)
    for i = 1:length(purviews_M)
        data.concept(concept2index(purviews_M{i},full_sys)).is_irreducible = 1;    
    end
end   