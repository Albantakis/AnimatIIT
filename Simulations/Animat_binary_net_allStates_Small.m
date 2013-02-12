clear all;
close all;

fprintf('\nRunning...\n\n');
tic;

%% Animat specificities
numSen = 2;
AnimatPath = '~/Documents/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_c1a3_24/trial1_';
%% options

%put options in array... TODO: this should be a struct and maybe also
%global

op_complex = 1;  % 0: only consider the whole system %1: find the complex
op_parallel = 0; % 1: parallel computing, 0: not
op_network = 4; % 1: random with self-connectivity, 2: random without self-connectivity, 3: modular, 4: logic gates, 0: some given connectivity matrix

op_TPM = 0; % 1: load TPM
op_ave = 1; % 0: use a specific current state 1: average over all possible current states
op_disp = 0; % 0: No figures, 1: only complex, 2: complex and whole system, 3: all figures
op_write = 1; % 1 --> write output to file
op_context = 0; % 0: conservative 1: progressive
op_empty = 1; % 0: excluding empty set in the past and the future 1: including empty set 
op_min = 1; % conservative only 0: phi is the sum of phi backward and phi forward (simulataneous partition)
                     % 1: phi is the minimum of phi_b and phi_f (separate partition)
op_console = 0; % 0: limited console output, 1: full console output
op_big_phi = 0; % 0 = big_phi is sum of small phi, 1 = big phi is volume best on EMD/best-packing, 2 = ave of H-difference b/w parent/child in hasse diagram
                     % 3 = look for pairwise distances less than 2*radius
                     % as marker for overlap
op_sum = 0; % 0 = compute big_phi_mip based on expanding parts into the space of the whole, 1 = just take whole minus sum of parts
op_normalize_big_phi = 1; % 0 = don't normalize big_phi, 1 = normalize big_phi but choose non-norm value, 2 = normalize but choose norm value
op_normalize_small_phi = 1; % 0 = don't normalize small_phi, 1 = normalize small_phi but choose non-norm value, 2 = normalize but choose norm value

%% inactive options, which are not used anymore
op_fb = 3; % 0: forward repertoire, 1: backward repertoire 2: both separately 3: both simultaneously
op_phi = 1; % two versions of small phi 0:Difference of entropy, 1:KL-divergence
op_whole = 0; % KLD is computed in 0: small system 1: whole system (previous version)

global grain, global noise;
grain = 100;
noise = 0.0; % 0 <= noise <= .5, this noise is applied to all nodes on their output essentially adding uncertaintity
global BRs, global FRs


fprintf('Noise Level = %f\n',noise);

options = [op_fb op_phi op_disp 1 1 op_context op_whole op_empty op_min op_console op_big_phi op_sum op_normalize_big_phi op_normalize_small_phi op_complex];

save options options

%% load tpm
%ToDo manually: used_nodes = Elements, current state, load correct mat file
%and correct edge list!
op_2TS = 0;    %do tpm after 2nd timestep
%for g = [261:10:441 442:461 471:10:901 902:911 921:10:1499]    %trial0_a24c13_s0
%for g = [547:10:727 737:807 817:10:1499]                       %trial1_a24c13_s0
%for g = [145:10:475 476:495 505:10:595 596:675 685:10:815 816:825 835:10:1499] %trial2_a24c13_s0
%for g = [173:10:213 214:222 223:10:313 314:323 333:10:443 444:463 473:10:663 664:683 683:10:1499]
%for g = [405:10:445 446:455 465:10:695 696:705 715:10:875 876:885 895 905 906:915 925:10:1499]
for g = 15392
 

Animat_gen = g;
global J

if op_TPM == 1
    used_nodes = [2     3     4     8    11    12    13    14    15    16];
    %used_nodes = [1     2     3     4     5     6     8    11    12    13    15    16];
    %used_nodes = [1,2,3,4,5,6,7,8,9,10,11,13,15,16];

    %load Animat1_tpm_short6;
    %load Animat157_8short_9to16Rest0;
    %load Animat157_6short_11to16Rest0;
    tpmfile = strcat('../Animat', int2str(Animat_gen), '_tpm');
    load(tpmfile)    
else
    [tpm, used_nodes] = LA_Animat_ReducedTpmSmall(Animat_gen, AnimatPath);
end
tpm(:,used_nodes < numSen+1) = 0.5; %Sensors might be switched on through mechanism but that is overwritten by environment --> doesn't do anything

J_tempfile = strcat(AnimatPath, int2str(Animat_gen), '_EdgeList.txt');
J_temp = load(J_tempfile);

TotNodes = max(max(J_temp))+1;
if ~isempty(TotNodes)

if op_write == 1
    Foldername = strcat('Animat_gen', int2str(Animat_gen));
    mkdir(Foldername)
end 
   
J = zeros(TotNodes);

for i = 1:TotNodes 
    out_ind = find(J_temp(:,1) == i-1); %minus one because nodes are labeled 0 to 15
    if ~isempty(out_ind)
        J(J_temp(out_ind,2)+1, i) = 1;
    end
    in_ind = find(J_temp(:,2) == i-1); %minus one because nodes are labeled 0 to 15
    if ~isempty(in_ind)
        J(i, J_temp(in_ind,1)+1) = 1;
    end
end

J = J(used_nodes, used_nodes);
J(used_nodes < numSen+1,:) = 0; %Sensors might have incoming connection, but because they actually don't do anything they shouldn't be taken into account

p = tpm;
p_x0 = p;
N = size(p,2);

BRs = cell(2^N,2^N); % backward repertoire
FRs = cell(2^N,2^N); % forward repertoire

% current state
if op_ave == 0
    %current_state = zeros(N,1); % all OFF
    %current_state = ones(N,1); % all ON
    %current_state = double(sum(tpm)>length(tpm)/2)';
    current_state = [0.5 0.5 0.5 0 0 0 1 1 1]'; 
    %current_state = [  0     0     0     0     1     0     1     1     1     1]';
    z_max = 1;
else
    [LifeStates, p_LifeStates] = LA_Animat_LifeStateListSmall(Animat_gen, used_nodes, AnimatPath, numSen);
    z_max = size(LifeStates, 1);
end
%% binary table
b_table = cell(2^N,N);
states = zeros(N,2^N);
for i=1: N
    for j=1: 2^i
        b_table{j,i} = trans2(j-1,i);
        if i== N
            states(:,j) = trans2(j-1,i);
        end
    end
end

%% parallel computing
isOpen = matlabpool('size');
if  isOpen == 0 && op_parallel > 0
%     s = ['matlabpool ' int2str(op_parallel)];
%     eval(s);
    matlabpool;
end

%% compute the big-phi
% if op_fb == 0
%     p = p'; % forward repertoire
% else
%     H_max = N;
% end

Big_phi_st = zeros(z_max,1);
Big_phi_MIP_st = zeros(z_max,1);

for z=1: z_max
    BRs = cell(2^N,2^N); % backward repertoire
    FRs = cell(2^N,2^N); % forward repertoire
    if op_write == 1
       diaryfile = strcat(Foldername, '/','Animat_gen', int2str(Animat_gen),'_state', int2str(z));
       diary(diaryfile); 
       diary on
    end    
    if op_ave == 0
        x1 = current_state;
    else
        x1 = LifeStates(z, :)';
        state = strcat('[', int2str(x1'),']')
    end
%     fprintf('x1=%s\n',mat2str(x1));
    
    % partial_prob_comp(partition, partition, state, prob_matrix, binary
    % table, op_fb
    check_prob = partial_prob_comp(1:N,1:N,x1,p,b_table,1); % last argument is op_fb = 1;
    state_check = sum(check_prob);
    if state_check == 0
        fprintf('This state cannot be realized!\n')
        Big_phi_st(z) = 0;
        Big_phi_MIP_st(z) = 0;
    else
        if op_complex == 0 % only consider whole system
            if op_fb == 2
                options(1) = 0; [Big_phi_f phi_f prob_cell_f] = big_phi_comp(1:N,x1,p,b_table,options);
                options(1) = 1; [Big_phi_b phi_b prob_cell_b] = big_phi_comp(1:N,x1,p,b_table,options);
                Big_phi = Big_phi_f + Big_phi_b;
            elseif op_fb == 3 % THIS IS THE ONLY ONE WE DO NOW? BOTH FORWARD AND BACKWARD SIMULTANEOUSLY
                M = 1:N;
                if op_context == 0
%                     [BRs FRs] = comp_pers(x1,p,b_table,options);
                    [Big_phi phi prob_cell MIPs M_IRR] = big_phi_comp_fb(M,x1,p,b_table,options);
                    % irreducible points
                    [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(prob_cell,phi,MIPs,M, 0,op_fb);
                    fprintf('\n')
                    fprintf('---------------------------------------------------------------------\n\n')
                    fprintf('Big_phi = %f\n', Big_phi);
                    fprintf('Sum of small_phis = %f\n',sum(phi));
                    fprintf('\nCore Concepts For Complex (Purview, MIP(past & future), Small phi):\n\n');
                    plot_REP(Big_phi, IRR_REP,IRR_phi,IRR_MIP, 1, M, options)
                else
                    [Big_phi phi prob_cell MIP prob_cell2] = big_phi_comp_fb(M,x1,p,b_table,options);
                end
            else
                [Big_phi phi prob_cell] = big_phi_comp(1:N,x1,p,b_table,options);
            end
            Big_phi_st(z) = Big_phi;
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THE CURRENT SETTINGS TAKE US HERE    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else % find the complex
            [Big_phi_MIP MIP Big_phi_M IRR_phi IRR_REP IRR_MIP M_IRR prob_M phi_M MIP_M] ...
                = big_phi_complex(x1,p,b_table,options);
            
            if op_fb == 2
                % subindex b means backward and f means forward
                IRR_phi_b = IRR_phi{1};
                IRR_phi_f = IRR_phi{2};
                IRR_REP_b = IRR_REP{1};
                IRR_REP_f = IRR_REP{2};
                M_IRR_b = M_IRR{1};
                M_IRR_f = M_IRR{2};
                % state dependent big phi and big phi MIP
                Big_phi_st(z) = sum(IRR_phi_b)+sum(IRR_phi_f);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THE CURRENT SETTINGS TAKE US HERE    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                Big_phi_st(z) = sum(IRR_phi);
            end
            Big_phi_MIP_st(z) = Big_phi_MIP;
        end
    end
    diary off
    % pause;
end

if op_ave == 1
    if op_fb == 0
        Big_phi_ave = sum(Big_phi_st)/z_max;
    else
        Big_phi_ave = sum(p_LifeStates .* Big_phi_st); %weighted ave/expected value
        Big_phi_MIP_ave = sum(p_LifeStates .* Big_phi_MIP_st);
    end
    if op_write == 1
       diaryfileEND = strcat(Foldername, '/','Animat_gen', int2str(Animat_gen),'_All.m');
       diary(diaryfileEND); 
       diary on
       cd(Foldername)
       save('LifeStates', 'LifeStates');
       save('p_LifeStates', 'p_LifeStates');
       save('Big_phi_st', 'Big_phi_st');
       save('Big_phi_MIP_st', 'Big_phi_MIP_st');
       cd ..
    end  
    fprintf('Big_phi_ave=%f\n',Big_phi_ave);
    fprintf('Big_phi_MIP_ave=%f\n',Big_phi_MIP_ave);
end

op_close = 0;
isOpen = matlabpool('size');
if isOpen > 0 && op_close == 1
    matlabpool close;
end

fprintf('\n');
if length(used_nodes) == 2
    used_nodes = used_nodes';   %Larissa: was occasionally necessary. not sure if it's still
end   
Element_mapping = [1:N; used_nodes-1]
fprintf('\n');

toc;

if op_write == 1
       diary off
end

%grPlot([], unique(J_temp, 'rows'), 'd');
end
end
toc;