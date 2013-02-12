clear all
Elem = 0:7;
plotflag = 3;
trialnum = [1 3 4 5 7];
numtrials = numel(trialnum);
totsteps = 30000-1;
step = 16;
range = 0:step:totsteps;
%startpoints = [696, 189, 349, 1692, 2356, 4372];
%cond = 'c1a3_24';
%cond = 'c24a35_36';
cond = 'c1b11a2b5_36'
path = strcat('~/Documents/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', cond, '/trial');
MaxFitness = 128;

BigPhi = zeros(numtrials, length(range));
BigPhiMip = zeros(numtrials,length(range));
BigPhiMip_B = zeros(numtrials,length(range));
BigPhiMip_max = zeros(numtrials,length(range));
diffStates = zeros(numtrials,length(range));
NumLifeStates = zeros(numtrials,length(range));
MeanNumConcepts = zeros(numtrials,length(range));
MeanSizeComplex = zeros(numtrials,length(range));
MeanHOConcepts = zeros(numtrials,length(range));
Fitness_level = zeros(numtrials,length(range));

ElemUsed = zeros(numtrials,numel(Elem),length(range));
ElemUsed_temp = zeros(1, numel(Elem));
FitPhiCorr = [];

cd Simulations
for t = 1:numtrials
    for i = 1:length(range)
        %------------- get Fitness from Animat files---------------------------
        docname = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_EdgeList.txt');
        docname2 = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');
        Fitness = load(docname2);
        Fitness_level(t,i) = Fitness(1);  
        EdgeList = load(docname);
        if ~isempty(EdgeList)
            Elements = unique(EdgeList); 
            ElemUsed_temp(Elements+1) = 5;
            Input = unique(EdgeList(:,1));
            Output = unique(EdgeList(:,2));
            ElemUsed_temp(setdiff(Elements, Output)+1) = 2;
            ElemUsed_temp(setdiff(Elements, Input)+1) = 1;
            ElemUsed(t,:,i) = ElemUsed_temp;
        end
      
        %------------- get rest from Phi calculation --------------------------
        % only go to folder here, because I need to access functions from the
        % main folder then, so I need to go out of the folder after loading
        cd(strcat(cond, '_trial', int2str(trialnum(t))))
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
        
        
        BigPhi(t,i) = sum(results.p_LifeStates .* Phi); %weighted ave/expected value
        BigPhiMip(t,i) = sum(results.p_LifeStates .* Phi_MIP);  
        BigPhiMip_B(t,i) = mean(Phi_MIP);  
        BigPhiMip_max(t,i) = max(Phi_MIP);  
        
        diffStates(t,i) = length(unique(Phi_MIP));
        NumLifeStates(t,i) = size(results.p_LifeStates, 1);
        MeanNumConcepts(t,i) = mean(Num_Concepts);
        MeanSizeComplex(t,i) = mean(subsys_numel);
        MeanHOConcepts(t,i) = mean(MConceptOrder);
        FitPhiCorr = [FitPhiCorr; [Fitness(1), BigPhiMip(t,i)]];
        else 
            cd ..
        end      
    end
end    
cd ..
%%
Fitness_Phi_Corr = corrcoef(mean(Fitness_level),mean(BigPhiMip))
MeanFitness = mean(Fitness_level,1);
MeanFitness = 100.*MeanFitness./MaxFitness;

if plotflag == 3
    figure
    subplot(3,1,1)
    hold on
    plot(range, MeanFitness)
    xlim([1, totsteps])
    
    subplot(3,1,2)
    hold on
    plot(range, mean(BigPhi), '-b')
    plot(range, mean(BigPhiMip_B), '-g')
    plot(range, mean(BigPhiMip_max), '-m')
    plot(range, mean(BigPhiMip), '-r')
    xlim([1, max(range)])

    subplot(3,1,3)
    hold on
    %plot(range, NumConcept, '.-b')
    plot(range, mean(NumLifeStates), '-b')
    plot(range, mean(diffStates), '-r')
    plot(range, mean(MeanNumConcepts), '-k')
    plot(range, mean(MeanSizeComplex), '-g')
    plot(range, mean(MeanHOConcepts), '-r')
    xlim([1, max(range)])
    
else
    figure(1)
    subplot(2,1,1)
    hold on
    plot(range, mean(BigPhi), '-b')
    plot(range, mean(BigPhiMip), '-r')
    xlim([1, max(range)])

    subplot(2,1,2)
    hold on
    %plot(range, NumConcept, '.-b')
    plot(range, mean(diffStates), '-r')
    xlim([1, max(range)])
end

% 
%  mBigPhi=[]; mBigPhiMip_B=[]; mBigPhiMip_max=[]; mBigPhiMip=[]; mNumLifeStates=[]; mdiffStates=[];
% for s = 0:50:totsteps-50 
%     mBigPhi = [mBigPhi mean(sum(BigPhi(:,s+1:s+50),2))];
%     mBigPhiMip_B = [mBigPhiMip_B mean(sum(BigPhiMip_B(:,s+1:s+50),2))];
%     mBigPhiMip_max = [mBigPhiMip_max mean(sum(BigPhiMip_max(:,s+1:s+50),2))];
%     mBigPhiMip = [mBigPhiMip mean(sum(BigPhiMip(:,s+1:s+50),2))];
%     mNumLifeStates = [mNumLifeStates mean(sum(NumLifeStates(:,s+1:s+50),2))];
%     mdiffStates = [mdiffStates mean(sum(diffStates(:,s+1:s+50),2))];
% end
% 
% if plotflag == 2
%     figure
%     subplot(3,1,1)
%     hold on
%     plot(mean(Fitness_level,1))
%     xlim([1, totsteps])
%     
%     subplot(3,1,2)
%     hold on
%     plot(25:50:totsteps, mBigPhi, '.-b')
%     plot(25:50:totsteps, mBigPhiMip_B, '.-g')
%     plot(25:50:totsteps, mBigPhiMip_max, '.-m')
%     plot(25:50:totsteps, mBigPhiMip, '.-r')
%     xlim([1, totsteps])
% 
%     subplot(3,1,3)
%     hold on
%     %plot(range, NumConcept, '.-b')
%     plot(25:50:totsteps, mNumLifeStates, '.-b')
%     plot(25:50:totsteps, mdiffStates, '.-r')
%     xlim([1, totsteps])
% end
% 
% 
if plotflag == 2
    figure
    subplot(2,1,1)
    hold on
    plot(FitPhiCorr(:,1), FitPhiCorr(:,2), '.k')
    
    [Fitind, Frow, Fcol] = unique(FitPhiCorr(:,1));
    for i= 1:length(Fitind)
        FPhi(i) = mean(FitPhiCorr(Fcol == i,2));
    end
    
    subplot(2,1,2)
    hold on
    plot(Fitind, FPhi)
end    
