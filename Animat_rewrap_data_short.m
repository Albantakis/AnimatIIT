function data = Animat_rewrap_data_short(Big_phi_M, phi_M, prob_M, M_cell, concept_MIP_M, purviews_M,...
                        Big_phi_MIP, MIP, Complex, M_i_max, BFCut)

                    
    %Only save main complex
data.Phi_MIP = Big_phi_MIP;
data.BFCut = BFCut;

    % store main complex of whole system as index
subsys = subsystem2index(Complex); 
data.main_complex = subsys; 
    
data.Phi = Big_phi_M(subsys);

num_concepts = 2^length(Complex) - 1;

% MIP for this subsystem
data.MIP1 = MIP(1);
data.MIP2 = MIP(2);


    % BUILD CONCEPTS!

% initialize size of concept struct array
data.concept(num_concepts).phi.min = 0;

for concept_i = 1:num_concepts

    % for now we set concept to reducible and correct after for-loop
    % (see below)
    data.concept(concept_i).is_irreducible = 0;
        
    if ~isempty(phi_M{subsys})
        % phi values for each concept
        data.concept(concept_i).phi.min = phi_M{subsys}(concept_i,1);
        data.concept(concept_i).phi.backwards = phi_M{subsys}(concept_i,2);
        data.concept(concept_i).phi.forwards = phi_M{subsys}(concept_i,3);

        %distribution
    %     data.concept(concept_i).distribution.backwards.whole = prob_M{subsys,1}{concept_i}{1};
    %     data.concept(concept_i).distribution.backwards.MIP = prob_M{subsys,2}{concept_i}{1};
    %     data.concept(concept_i).distribution.forwards.whole = prob_M{subsys,1}{concept_i}{2};
    %     data.concept(concept_i).distribution.forwards.MIP = prob_M{subsys,2}{concept_i}{2};
    % 
        % denominator
        data.concept(concept_i).denominator.backwards = ...
                                                            sort([concept_MIP_M{subsys}{concept_i}{:,1,1}]);
        data.concept(concept_i).denominator.forwards = ...
                                                            sort([concept_MIP_M{subsys}{concept_i}{:,1,2}]);

        % MIPs
    %     data.concept(concept_i).MIP.backwards.numerator1 = sort([concept_MIP_M{subsys}{concept_i}{1,2,1}]);
    %     data.concept(concept_i).MIP.backwards.numerator2 = sort([concept_MIP_M{subsys}{concept_i}{2,2,1}]);
    %     data.concept(concept_i).MIP.backwards.denominator1 = sort([concept_MIP_M{subsys}{concept_i}{1,1,1}]);
    %     data.concept(concept_i).MIP.backwards.denominator2 = sort([concept_MIP_M{subsys}{concept_i}{2,1,1}]);
    %     data.concept(concept_i).MIP.forwards.numerator1 = sort([concept_MIP_M{subsys}{concept_i}{1,2,2}]);
    %     data.concept(concept_i).MIP.forwards.numerator2 = sort([concept_MIP_M{subsys}{concept_i}{2,2,2}]);
    %     data.concept(concept_i).MIP.forwards.denominator1 = sort([concept_MIP_M{subsys}{concept_i}{1,1,2}]);
    %     data.concept(concept_i).MIP.forwards.denominator2 = sort([concept_MIP_M{subsys}{concept_i}{2,1,2}]);
    end
end

% label irreducible concepts
if ~isempty(purviews_M{subsys})

    for i = 1:length(purviews_M{subsys})

        data.concept(concept2index(purviews_M{subsys}{i},index2subsystem(subsys))).is_irreducible = 1;    

    end
end
end    

    
    
    
