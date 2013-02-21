function data = rewrap_data(Big_phi_M, phi_M, prob_M, M_cell, concept_MIP_M, purviews_M,...
                        Big_phi_MIP, MIP, Complex, M_i_max,  Big_phi_MIP_M, complex_MIP_M, Big_phi_MIP_all_M, complex_MIP_M_all)
                    
                    

N_subsys = length(Big_phi_M);
data.subsystem(N_subsys).Phi = 0;

% store main complex of whole system as index
data.main_complex = subsystem2index(Complex);

% Phi values for each subset, "Big Phi"
temp_cell = mat2cell(Big_phi_M,ones(1,N_subsys),1);
[data.subsystem.Phi] = temp_cell{:};

% Phi_MIP values for each subset, "Big Phi MIP"
temp_cell = mat2cell(Big_phi_MIP_M,ones(1,N_subsys),1);
[data.subsystem.Phi_MIP] = temp_cell{:};


    

for subsys = 1:N_subsys

    if ~isempty(purviews_M{subsys})

        subsystem = M_cell{subsys};
        num_concepts = 2^length(subsystem) - 1;

        % MIP for this subsystem
        data.subsystem(subsys).MIP1 = complex_MIP_M{subsys};
        data.subsystem(subsys).MIP2 = pick_rest(subsystem,complex_MIP_M{subsys});

        % store Big_Phi_Partitioned values for other partitions of this
        % sybsystem as well as a way to see what the partitions were... this is
        % the weakest link as far as storage... tackle another day...
        data.subsystem(subsys).partition_values = Big_phi_MIP_all_M{subsys};
        data.subsystem(subsys).partitions = complex_MIP_M_all{subsys};


        % BUILD CONCEPTS!

        % initialize size of concept struct array
        data.subsystem(subsys).concept(num_concepts).phi.min = 0;

        for concept_i = 1:num_concepts

            % for now we set concept to reducible and correct after for-loop
            % (see below)
            data.subsystem(subsys).concept(concept_i).is_irreducible = 0;

            % phi values for each concept
            data.subsystem(subsys).concept(concept_i).phi.min = phi_M{subsys}(concept_i,1);
            data.subsystem(subsys).concept(concept_i).phi.backwards = phi_M{subsys}(concept_i,2);
            data.subsystem(subsys).concept(concept_i).phi.forwards = phi_M{subsys}(concept_i,3);

            %distribution
            data.subsystem(subsys).concept(concept_i).distribution.backwards.whole = prob_M{subsys,1}{concept_i}{1};
            data.subsystem(subsys).concept(concept_i).distribution.backwards.MIP = prob_M{subsys,2}{concept_i}{1};
            data.subsystem(subsys).concept(concept_i).distribution.forwards.whole = prob_M{subsys,1}{concept_i}{2};
            data.subsystem(subsys).concept(concept_i).distribution.forwards.MIP = prob_M{subsys,2}{concept_i}{2};

            % denominator
            data.subsystem(subsys).concept(concept_i).denominator.backwards = ...
                                                                sort([concept_MIP_M{subsys}{concept_i}{:,1,1}]);
            data.subsystem(subsys).concept(concept_i).denominator.forwards = ...
                                                                sort([concept_MIP_M{subsys}{concept_i}{:,1,2}]);

            % MIPs
            data.subsystem(subsys).concept(concept_i).MIP.backwards.numerator1 = sort([concept_MIP_M{subsys}{concept_i}{1,2,1}]);
            data.subsystem(subsys).concept(concept_i).MIP.backwards.numerator2 = sort([concept_MIP_M{subsys}{concept_i}{2,2,1}]);
            data.subsystem(subsys).concept(concept_i).MIP.backwards.denominator1 = sort([concept_MIP_M{subsys}{concept_i}{1,1,1}]);
            data.subsystem(subsys).concept(concept_i).MIP.backwards.denominator2 = sort([concept_MIP_M{subsys}{concept_i}{2,1,1}]);
            data.subsystem(subsys).concept(concept_i).MIP.forwards.numerator1 = sort([concept_MIP_M{subsys}{concept_i}{1,2,2}]);
            data.subsystem(subsys).concept(concept_i).MIP.forwards.numerator2 = sort([concept_MIP_M{subsys}{concept_i}{2,2,2}]);
            data.subsystem(subsys).concept(concept_i).MIP.forwards.denominator1 = sort([concept_MIP_M{subsys}{concept_i}{1,1,2}]);
            data.subsystem(subsys).concept(concept_i).MIP.forwards.denominator2 = sort([concept_MIP_M{subsys}{concept_i}{2,1,2}]);

        end

        % label irreducible concepts
        if ~isempty(purviews_M{subsys})

            for i = 1:length(purviews_M{subsys})

                data.subsystem(subsys).concept(concept2index(purviews_M{subsys}{i},index2subsystem(subsys))).is_irreducible = 1;    

            end
        end
    end   
end
    
    
    
