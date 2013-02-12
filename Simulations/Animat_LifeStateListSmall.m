function [LifeStates, p_LifeStates, InStates, p_InStates] = Animat_LifeStateListSmall(gen, Elements, APath, numSen)
%LA 10/05/2012: Change: Set Sensors to .5 before finding unique states!

if nargin < 3 
    APath = '~/Documents/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work/trial0_';
end
if nargin < 1
    gen = 17;
    %Elements = [2     3     4     8    11    12    13    14    15    16];
end
%Creates List of states the animat gets into during its life plus vector of
%state probability
lifestatefile = strcat(APath, int2str(gen), '_LifetimeLogicTable.txt');
Life_tpm = importdata(lifestatefile,',', 1);

%% Environment (Input) distribution
% calculate probability distribution of input states
Input_vec = Life_tpm.data(:,1:numSen);
[InStates, inR, inMap] = unique(Input_vec, 'rows', 'First');
% all possible input states should appear in life states because of
% initialization of program
InStates_distr = hist(inMap, 1:size(InStates,1))';
p_InStates = InStates_distr./size(Input_vec,1);

%second order LifeStates correlations
Input_O2 = [Input_vec(1:end-1,:), Input_vec(2:end,:)];
[InStatesO2, inR, inMapO2] = unique(Input_O2, 'rows', 'First');
InStates_distrO2 = hist(inMapO2, 1:size(InStatesO2,1))';
p_InStates = InStates_distrO2./size(Input_O2,1);

%% Internal state distribution
Life_tpm = Life_tpm.data(:,10:17);
Life_tpm = Life_tpm(:,Elements); %reduced logic table
Life_tpm(:,Elements < numSen+1) = 0; %LA: For the moment but check later, maybe they are already anyways. Just to be safe
[LifeStates, redR, redMap] = unique(Life_tpm, 'rows', 'First');

LifeStates_distr = hist(redMap, 1:size(LifeStates,1))'; %Larissa Check if number in hist fits with LifeState Ordering!!!
p_LifeStates = LifeStates_distr./size(Life_tpm,1);

%B = redLogic(sortind,:); This is the correct ordering of input states
end