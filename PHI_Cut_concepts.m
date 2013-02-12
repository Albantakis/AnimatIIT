function [CutDistr] = PHI_Cut_concepts(M,MIP,BFCut,M_IRR_M,prob_M, phi_M,concept_MIP_M, network)
%%
% Recompute distributrions of PHI cut
%%
op_strongconn = network.options(10);
whole_i = subsystem2index(M);     % Full system index

phi_whole = phi_M{whole_i}(:,1)';
concept_numind = find(phi_whole ~= 0);
phi_w_concepts = phi_whole(phi_whole ~= 0);
IRR_whole = M_IRR_M{whole_i};

M1 = MIP;
M2 = pick_rest(M,M1);

M1_i = subsystem2index(M1);
M2_i = subsystem2index(M2);  

if BFCut == 1 %M1 <- M2 is cut
    BRcut_dist = cell(length(phi_w_concepts), 2, 2); %dim1: per concept, dim2: past/future, dim2: whole/cut
    for k = 1:length(phi_w_concepts)
        IRR_w = IRR_whole{k};
        if all(ismember(IRR_w,M1)) && op_strongconn ~= 0 %~isempty(prob_M{M1_i,1})
            % for M1 <- M2 cut take BR of M1 and FR from M
            indm = concept2index(IRR_w,M1);
            cutpdist = expand_prob(prob_M{M1_i,1}{indm}{1},M,M1);
            BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} cutpdist};
            BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}}; 
        elseif all(ismember(IRR_w,M2))  && op_strongconn ~= 0 % ~isempty(prob_M{M2_i,1})
            indm = concept2index(IRR_w,M2);
            %Larissa: distributions that are identical anyways

            %compute the max ent forward dist (or the marginal forward) of M1
            forward_max_ent_M1 = comp_pers_cpt(network.nodes,[],M1,[],'forward');
            cutfdist = expand_prob_general(prob_M{M2_i,1}{indm}{2},M,M2,forward_max_ent_M1(:));
            BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}}; %back is the same
            BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} cutfdist}; %future might have changed
        else % if numerator has elements from both sides
            denom_p = sort([concept_MIP_M{whole_i}{concept_numind(k)}{:,1,1}]);    %Larissa: The sort may be important 
            denom_f = sort([concept_MIP_M{whole_i}{concept_numind(k)}{:,1,2}]);
            % for BRcut (M1 <- M2 is cut) M1M2/[M1]p[M2]f is still intact
            BRcut_pdist = []; BRcut_fdist = []; 
            if all(ismember(denom_p,M1))
                if ~all(ismember(denom_f,M2))
                    % check the new forward phi_mip M1M2/[M1M2]f for all
                    % possible denominators
                    [phi_BRcut_FR, BRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','BRcut');
                end    
            else
                if all(ismember(denom_f,M2))
                    % check the new backward phi_mip M1M2/[M1M2]p for all
                    % possible denominators
                    [phi_BRcut_BR, BRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','BRcut');
                else
                    % check the new back and forward phi_mip M1M2/[M1M2]p[M1M2]f for all
                    % possible denominators
                    [phi_BRcut_BR, BRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','BRcut');
                    [phi_BRcut_FR, BRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','BRcut');
                end    
            end

            if ~isempty(BRcut_pdist)
                BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} BRcut_pdist};
            else    
                BRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}};
            end
            if ~isempty(BRcut_fdist)
                BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} BRcut_fdist};
            else
                BRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}};
            end
        end % if 
    end %for k    
    CutDistr = BRcut_dist;
else %M1 -> M2 is cut
    FRcut_dist = cell(length(phi_w_concepts), 2, 2);
    for k = 1:length(phi_w_concepts)
        IRR_w = IRR_whole{k};
        if all(ismember(IRR_w,M1)) 
            % for M1 <- M2 cut take BR of M1 and FR from M
            indm = concept2index(IRR_w,M1);
            %Larissa: distributions that are identical anyways are empty
            %compute the max ent forward dist (or the marginal forward) of M2
            forward_max_ent_M2 = comp_pers_cpt(network.nodes,[],M2,[],'forward');
            cutfdist = expand_prob_general(prob_M{M1_i,1}{indm}{2},M,M1,forward_max_ent_M2(:)); 
            FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}}; 
            FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} cutfdist};
        elseif all(ismember(IRR_w,M2))
            indm = concept2index(IRR_w,M2);
            %Larissa: distributions that are identical anyways
            %are empty for option 1
            cutpdist = expand_prob(prob_M{M2_i,1}{indm}{1},M,M2); 
            FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} cutpdist}; %back might have changed, future is the same
            FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}}; %back might have changed, future is the same
        else % if numerator has elements from both sides
            denom_p = sort([concept_MIP_M{whole_i}{concept_numind(k)}{:,1,1}]);    %Larissa: The sort may be important 
            denom_f = sort([concept_MIP_M{whole_i}{concept_numind(k)}{:,1,2}]);
            % for BRcut (M1 <- M2 is cut) M1M2/[M1]p[M2]f is still intact
            FRcut_pdist = []; FRcut_fdist = []; %Only needed for op_big_phi = 6 or 7

            % for FRcut (M1 -> M2 is cut) M1M2/[M2]p[M1]f is still intact
            if all(ismember(denom_p,M2))
                if ~all(ismember(denom_f,M1))
                    % check the new forward phi_mip M1M2/[M1M2]f for all
                    % possible denominators
                    [phi_FRcut_FR, FRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','FRcut');
                end    
            else
                if all(ismember(denom_f,M1))
                    % check the new backward phi_mip M1M2/[M1M2]p for all
                    % possible denominators
                    [phi_FRcut_BR, FRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','FRcut');
                else
                    % check the new back and forward phi_mip M1M2/[M1M2]p[M1M2]f for all
                    % possible denominators
                    [phi_FRcut_BR, FRcut_pdist, denom_pnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'backward','FRcut');
                    [phi_FRcut_FR, FRcut_fdist, denom_fnew, network] = phi_comp_ex_unidir(M,M1,M2,IRR_w,network.current_state,network,'forward','FRcut');
                end    
            end

            if ~isempty(FRcut_pdist)
                FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} FRcut_pdist};
            else
                FRcut_dist(k,1,:) = {prob_M{whole_i,1}{concept_numind(k)}{1} prob_M{whole_i,1}{concept_numind(k)}{1}};
            end
            if ~isempty(FRcut_fdist)
                FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} FRcut_fdist};
            else
                FRcut_dist(k,2,:) = {prob_M{whole_i,1}{concept_numind(k)}{2} prob_M{whole_i,1}{concept_numind(k)}{2}};
            end  
        end % if 
    end %for k 
    CutDistr = FRcut_dist;
end

%%
figure
L = length(CutDistr{1,1,1});
rows = length(phi_w_concepts);
states = convert(L);
%title('[1]->[23]');
for k = 1:rows
    subplot(rows,2,2*k-1)
        bar([CutDistr{k,1,1} CutDistr{k,1,2}])
        axis([0.5 L+0.5 0 1])
        set(gca,'XTick',1:L)
        set(gca,'YTickLabel',[0 .5 1])
        title(strcat(mat2str(IRR_whole{k}),'p'));
        if k == rows
            set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
            rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
        else
            set(gca,'XTickLabel',[]) 
        end
        
    subplot(rows,2,2*k)
        bar([CutDistr{k,2,1} CutDistr{k,2,2}])
        axis([0.5 L+0.5 0 1])
        set(gca,'XTick',1:L)
        set(gca,'YTickLabel',[0 .5 1])
        title(strcat(mat2str(IRR_whole{k}),'f'));

        if k == rows
            set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
            rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
        else
            set(gca,'XTickLabel',[]) 
        end

end    
end

%%
function states = convert(N)
states = zeros(N,log2(N));
for i=1: N
    sigma = trans2(i-1,log2(N));
    states(i,:) = sigma';
end
end