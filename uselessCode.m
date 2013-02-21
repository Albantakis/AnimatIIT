%% Find strongly connected system parts
if op_complex == 1
    % Check PHI for potential complexes and self-loops
    % This is a first sweep, later there is another check if the
    % subsystem tested is strongly connected
    FullSys = 1:num_nodes;
    SelfLoopElements = FullSys(diag(in_connect_mat)==1);
    if ~isempty(SelfLoopElements)
        for i = 1:length(SelfLoopElements)
            indexSelfLoop(i) = subsystem2index(SelfLoopElements(i));
        end
    end    

    J_sparse = sparse(in_connect_mat);

    [X,whichPosComp] = graphconncomp(J_sparse);

    PosComplexes = unique(whichPosComp);
    countComplex = hist(whichPosComp, PosComplexes);
    indexToRepeatedValue = (countComplex~=1);
    PosComplexes = PosComplexes(indexToRepeatedValue);
    indexComplex = [];
    if ~isempty(PosComplexes)
        for i = 1:length(PosComplexes)
            TryComplex = FullSys(whichPosComp == PosComplexes(i));
            indexComplex(i) = subsystem2index(TryComplex);
        end
    end

    indexComplex = union(indexComplex, indexSelfLoop);
    indexComplex = sort(indexComplex);     % Only non-overlapping (!) subsystems that have to be searched! But: Subsystems of these can also exist
end