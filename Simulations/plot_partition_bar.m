% function plot_partition_bar(Big_phi, REP_cell,phi,MIP_cell, M, plot_panel)
function plot_partition_bar(concepts, nWholeConcepts, whole_phis, parts_phis, selected_indices, plot_panel,...
                            whole_purviews, parts_purviews, Big_phi_MIP)


nPartConcepts = length(parts_phis); 


% if op_figures
%     figure(20+fig_st);
%     plot_phi(phi,MIP_cell, M)
%     title('\phi  For Each Purview','FontSize',20);
%     xlabel('Purview','FontSize',16)
%     ylabel('\phi','FontSize',20)
% end

if isempty(selected_indices)
    whole_concepts = concepts(1:nWholeConcepts,:);
    parts_concepts = concepts(nWholeConcepts+1:end,:);
else    
    whole_concepts = concepts(selected_indices(selected_indices <= nWholeConcepts),:);
    parts_concepts = concepts(selected_indices(selected_indices > nWholeConcepts)-nWholeConcepts,:);
end

nWholeConcepts = size(whole_concepts,1);
nPartConcepts = size(parts_concepts,1);
N = max(nWholeConcepts,nPartConcepts);        

r = N;
c = 2;



% pw = 300; % panel width including margin in pixels
% ph = 50; % panel height including margin in pixels
% mb = 20; % bottom margin
% mt = 20; % top margin
% mh = 25; % margin height
% mw = 10; % margin width
% % panel_size = [0, 0, pw, ph]'; % panel size including margin
% panel_size_w = [50, 0, pw-mw*2, ph-mh]'; % panel size without margin
% fig_size = [0, 0,  pw*c + 100, ph*r + mb+mt + 50]'; % size of figure window
% pos_fig = [100,100,0,0]' + fig_size; % position of figure

if r >= 4
    scaling = r/4;
else
    scaling = 1;
end


if scaling > 1
      
    panel_pos = get(plot_panel,'Position');
    
    panel_height = panel_pos(4)*scaling;
    panel_pos(4) = panel_height;
    panel_pos(2) = 1 - panel_pos(4);
    
    set(plot_panel,'Position',panel_pos);
    
    left_marg_col1 = .03;
    left_marg_col2 = .57;
    top_start = 1 - .45/panel_height;
    vert_marg = .1/panel_height;
    plot_h = .075/panel_height;
    plot_w = .4;
    
%     setappdata(plot_panel,'limit',top_start-(vert_marg + plot_h)*(r-1)+.03)
    
else
    
    panel_height = 1;
    left_marg_col1 = .03;
    left_marg_col2 = .57;
    top_start = .55;
    vert_marg = .1;
    plot_h = .075;
    plot_w = .4;
    panel_pos = get(plot_panel,'Position');
    setappdata(plot_panel,'limit',panel_pos(2))
end

pos_vec = zeros(4,r,c);

for i=1: r
    for j=1:c
        if (j == 1)
            pos_vec(:,i,j) = [left_marg_col1 top_start-(vert_marg + plot_h)*(i-1) plot_w plot_h];
        else
            pos_vec(:,i,j) = [left_marg_col2 top_start-(vert_marg + plot_h)*(i-1) plot_w plot_h];
        end
    end
end

sy = ['Big Phi (partition) = ',num2str(Big_phi_MIP)];
    
    
for i = 1:N
        
        pos_whole = pos_vec(:,i,1)';
        pos_part = pos_vec(:,i,2)';
            
        if i <= nWholeConcepts
        
            prob_whole = whole_concepts(i,:);
            phi_whole = whole_phis(i);
            
        
            whole_title = [mod_mat2str(whole_purviews{i}) ',  \phi=' num2str(phi_whole)];

            if i == nWholeConcepts
                labelON = 1;
            else
                labelON = 0;
            end
            
            plot_BRFR(prob_whole,pos_whole,whole_title,labelON,plot_panel,'b');
            
        end
        
        if i <= nPartConcepts
            
            prob_part = parts_concepts(i,:);
            phi_part = parts_phis(i);
            
        
            part_title = [mod_mat2str(parts_purviews{i}) ',  \phi=' num2str(phi_part)];

            if i == nPartConcepts
                labelON = 1;
            else
                labelON = 0;
            end
            
            plot_BRFR(prob_part,pos_part,part_title,labelON,plot_panel,'r')
            
        end

        if i == N

            uicontrol('Style', 'text',...
            'String', {'Whole and Partitioned Concepts',sy},... %replace something with the text you want
            'Units','normalized',...
            'FontSize',16,...
            'BackgroundColor','w',...
            'Position', [0.33 1-.2/panel_height 0.33 0.12/panel_height],...
            'Parent',plot_panel,...
            'Clipping','on');

        end


end

end

                


function [] = plot_BRFR(concept,pos,plot_title,labelON,plot_panel,color)


subplot('Position',pos,'Parent',plot_panel)
h = bar(0:length(concept)-1,concept,'hist');
set(h,'FaceColor',color)
set(gca,'XTick',0:length(concept)-1)
set(gca,'YTickLabel',[0 .5 1])

states = convert(length(concept));
axis([-0.5 length(concept)-0.5 0 1.0])
title(plot_title)
if labelON== 1
    set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
    rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
%     xlab = xlabel('State (Node Order: [1...N])','FontSize',14,'Units','pixels');
%     disp(get(xlab,'Position'))
%     set(xlab,'Position',get(xlab,'Position') - [0 100 0])
% disp(get(xlab,'Position'))
else
    set(gca,'XTickLabel',[]) 
end


end



function states = convert(N)
states = zeros(N,log2(N));
for i=1: N
    sigma = trans2(i-1,log2(N));
    states(i,:) = sigma';
end
end
