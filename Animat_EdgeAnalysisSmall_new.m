range = 0:16:30000-1;
Elem = 0:7;
trial = 5;
cond = 'c24a35_36';

path = strcat('~/Documents/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_', cond, '/trial');
%path = '/Volumes/Macintosh' HD 2'/Postdoc/Projects/Hintze_Animat/temporalSpatialIntegrationLite/work/trial';
Fitness_level = zeros(1,length(range));
ElemUsed = zeros(numel(Elem),length(range));
for i = 1:length(range)
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
end

FitInc = find(diff(Fitness_level))+1;
figure;
subplot(2,1,1)
hold on
imagesc(ElemUsed);
for j = FitInc
    line([j j],[0 9], 'color', [1 0 0])
end  
axis([1 length(range) 1 8])
subplot(2,1,2)
hold on
plot(range, Fitness_level)
for j = FitInc
    line([range(j) range(j)],[min(Fitness_level) max(Fitness_level)], 'color', [1 0 0])
end  
xlim([1, max(range)])


%% Find the most common states
gen = 17;
docname3 = strcat('~/Documents/Arend_XCodeAnimat/representationCode/work/generation_', int2str(gen), '_LifetimeLogicTable.txt');
LifeLogic = importdata(docname3, ',', 1);
LifeLogic = LifeLogic.data;
[temp, a, b] = unique(LifeLogic(:,[18:33]), 'rows');
for i = 1: size(a)
    state_Freq(i) = nnz(b == i);
end    
[maxSt, maxStNb] = max(state_Freq);
state1 = LifeLogic(a(maxStNb),:); 

ind7 = find(LifeLogic(:,25)==1); % node 7 on 
LifeLogic_7 = LifeLogic(ind7, :);
[LL7, a7, b7] = unique(LifeLogic_7, 'rows');
%% Plot Life Logic Table
figure
    subplot(2,1,1)
    imagesc(LifeLogic(:, [1:16])')
    subplot(2,1,2)
    imagesc(LifeLogic(:, [18:33])')
    colormap('gray')