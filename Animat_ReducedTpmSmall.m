function [tpm, used_nodes, RedStates] = Animat_ReducedTpmSmall(gen, APath, numSen, numMot, J)

% NODES WITHOUT INPUTS OR OUTPUTS (except sensors) are taken out of the tpm
% NODES WITH ONLY OUTPUTS are set to 0 in animat program -> they shouldn't be
% marginalized -> take them out of TPM
% MOTORS are set to 0 after each step in the animat program
% STRATEGY: leave motors in TPM but always use 00 state for them
% In principle, IIT with removals should take care of nodes without in, or
% outputs, but previous exlusion doesn't hurt and might help speed

if nargin < 5 %like in previous versions
    J = [];
    docname2 = strcat(APath, int2str(gen), '_EdgeList.txt');
    EdgeList = load(docname2);
    used_nodes = (unique(EdgeList)+1)'; 
end

numNodes = size(J,1);

% Load Logic Table
docname = strcat(APath, int2str(gen), '_FullLogicTable.txt');
FullLogic = importdata(docname,',', 1);
FullLogic = FullLogic.data(:,[1:numNodes,(numNodes+2):(2*numNodes+1)]);

if ~isempty(J)
    %OnlyOut also finds nodes that have no input and output
    SInputs = sum(J, 2);
    OnlyOut = find(SInputs == 0);
    OnlyOut = OnlyOut(OnlyOut > numSen);
    %also exclude only input nodes that are not motors (for the moment at
    %least)
    SInputs = sum(J,1);
    OnlyIn = find(SInputs == 0);
    OnlyIn = OnlyIn(OnlyIn <= (numNodes-numMot));
    OnlyInOut = union(OnlyIn, OnlyOut);
    used_nodes = setdiff(1:numNodes,OnlyInOut);
end    

%Set motor-TPM to motor = all zeros
%indMot0 = find(sum(FullLogic(:,(numNodes-numMot+1):numNodes),2) == 0);
indMot0 = 1:2^(numNodes-numMot); %Motors are always last in TPM
tpm = repmat(FullLogic(indMot0,numNodes+(1:numNodes)),2^numMot,1);
FullLogic = [FullLogic(:,1:numNodes), tpm];

%Exclude nodes with only in or outputs (only Outputs can be dangerous)
excluded_nodes = setdiff(1:numNodes,used_nodes);
%take TPM for all excluded elements = 0. 0 important for OnlyOut Nodes, for
%not used nodes it wouldn't matter if 0 or 1 is used 
ind0 = find(sum(FullLogic(:,excluded_nodes),2) == 0);
RedStates = FullLogic(ind0,used_nodes); 
tpm = tpm(ind0,used_nodes);    
end