clear all
Elem = 0:7;
plotflag = 2;
range = 0:16:30000-1;
%cond = 'c1a3_24';
%cond = 'c24a35_36';
cond = 'c1b11a2b5_36';
trial = 9;

path = strcat('~/Documents/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', cond, '/trial');

Fitness_level = zeros(1,length(range));
ElemUsed = zeros(numel(Elem),length(range));

cd Simulations
for i = 1:length(range)
    %------------- get Fitness from Animat files---------------------------
    docname = strcat(path, int2str(trial), '_', int2str(range(i)), '_EdgeList.txt');
    docname2 = strcat(path, int2str(trial), '_', int2str(range(i)), '_KOdata.txt');
    Fitness = load(docname2);
    Fitness_level(i) = Fitness(1);  
    EdgeList = load(docname);
    if ~isempty(EdgeList)
        Elements = unique(EdgeList); 
        ElemUsed(Elements+1,i) = 5;
        Input = unique(EdgeList(:,1));
        Output = unique(EdgeList(:,2));
        ElemUsed(setdiff(Elements, Output)+1,i) = 2;
        ElemUsed(setdiff(Elements, Input)+1,i) = 1;
    end    
    
    %------------- get rest from Phi calculation --------------------------
    % only go to folder here, because I need to access functions from the
    % main folder then, so I need to go out of the folder after loading
    cd(strcat(cond, '_trial', int2str(trial)))
    Animat_gen = range(i);
    FilenameA = strcat('Animat_gen', int2str(Animat_gen),'_ShortResults.mat');
    FilenameB = strcat('Animat_gen', int2str(Animat_gen),'_Results.mat');
    
    if exist(FilenameA,'file') == 2 
        load(FilenameA)
        cd ..
        LSindex = results.LifeState_index;
        % initialize 
        subsys_numel = zeros(numel(LSindex),1);
        Phi = zeros(numel(LSindex),1);
        Phi_MIP = zeros(numel(LSindex),1);
        Num_Concepts = zeros(numel(LSindex),1);
        MConceptOrder = zeros(numel(LSindex),1);
        for j = 1:numel(LSindex)
            subsys = results.Complex{LSindex(j)};
            subsys_numel(j) = numel(subsys);
            subsys_ind = subsystem2index(subsys);
            sys_results = results.state(LSindex(j));
            Phi(j) = sys_results.Phi;
            Phi_MIP(j) = sys_results.Phi_MIP;
            IrrConcepts = [sys_results.concept(:).is_irreducible];
            Num_Concepts(j) = nnz(IrrConcepts);
            % higher level concepts
            numeratorindex = find(IrrConcepts);
            if Num_Concepts(j) > 0
                ConceptOrder = zeros(size(numeratorindex));
                for k = 1:length(numeratorindex)
                     ConceptOrder(k) = numel(index2concept(numeratorindex(k), subsys)) > 1;
                end
                MConceptOrder(j) = sum(ConceptOrder);
            end
        end
        BigPhi(i) = sum(results.p_LifeStates .* Phi); %weighted ave/expected value
        BigPhiMip(i) = sum(results.p_LifeStates .* Phi_MIP);  
        BigPhiMip_B(i) = mean(Phi_MIP);  
        BigPhiMip_max(i) = max(Phi_MIP);  
        
        diffStates(i) = length(unique(Phi_MIP));
        NumLifeStates(i) = size(results.p_LifeStates, 1);
        MeanNumConcepts(i) = mean(Num_Concepts);
        MeanSizeComplex(i) = mean(subsys_numel);
        MeanHOConcepts(i) = mean(MConceptOrder);
    elseif exist(FilenameB, 'file') == 2    
        load(FilenameB)
        cd ..
        LSindex = results.LifeState_index;
        % initialize 
        subsys_numel = zeros(numel(LSindex),1);
        Phi = zeros(numel(LSindex),1);
        Phi_MIP = zeros(numel(LSindex),1);
        Num_Concepts = zeros(numel(LSindex),1);
        MConceptOrder = zeros(numel(LSindex),1);
        for j = 1:numel(LSindex)
            subsys = results.Complex{LSindex(j)};
            subsys_numel(j) = numel(subsys);
            subsys_ind = subsystem2index(subsys);
            sys_results = results.state(LSindex(j)).subsystem(subsys_ind);
            Phi(j) = sys_results.Phi;
            Phi_MIP(j) = sys_results.Phi_MIP;
            IrrConcepts = [sys_results.concept(:).is_irreducible];
            Num_Concepts(j) = nnz(IrrConcepts);
            % higher level concepts
            numeratorindex = find(IrrConcepts);
            if Num_Concepts(j) > 0
                ConceptOrder = zeros(size(numeratorindex));
                for k = 1:length(numeratorindex)
                     ConceptOrder(k) = numel(index2concept(numeratorindex(k), subsys)) > 1;
                end
                MConceptOrder(j) = sum(ConceptOrder);
            end
        end
        BigPhi(i) = sum(results.p_LifeStates .* Phi); %weighted ave/expected value
        BigPhiMip(i) = sum(results.p_LifeStates .* Phi_MIP);  
        BigPhiMip_B(i) = mean(Phi_MIP);  
        BigPhiMip_max(i) = max(Phi_MIP);  
        
        diffStates(i) = length(unique(Phi_MIP));
        NumLifeStates(i) = size(results.p_LifeStates, 1);
        MeanNumConcepts(i) = mean(Num_Concepts);
        MeanSizeComplex(i) = mean(subsys_numel);
        MeanHOConcepts(i) = mean(MConceptOrder);
    else
        cd ..
    end    
end    
cd ..
%%

if plotflag == 2
    FitInc = find(diff(Fitness_level))+1;
    figure
    subplot(3,1,1)
    hold on
    plot(range, Fitness_level)
    for j = FitInc
        line([range(j) range(j)],[min(Fitness_level) max(Fitness_level)], 'color', [1 0 0])
    end  
    xlim([1, max(range)])
    
    %range = 0:16:1567;
    subplot(3,1,2)
    hold on
    plot(range, BigPhi, '-b')
    plot(range, BigPhiMip_B, '-g')
    plot(range, BigPhiMip_max, '-m')
    plot(range, BigPhiMip, '-r')
    xlim([1, max(range)])

    subplot(3,1,3)
    hold on
    %plot(range, NumConcept, '.-b')
    plot(range, NumLifeStates, '-b')
    plot(range, diffStates, '-r')
    plot(range, MeanNumConcepts, '-k')
    plot(range, MeanSizeComplex, '-g')
    plot(range, MeanHOConcepts, '-r')
    xlim([1, max(range)])
    
else
    figure(1)
    subplot(2,1,1)
    hold on
    plot(range, BigPhi, '-b')
    plot(range, BigPhiMip, '-r')
    xlim([1, max(range)])

    subplot(2,1,2)
    hold on
    %plot(range, NumConcept, '.-b')
    plot(range, diffStates, '-r')
    xlim([1, max(range)])
end

