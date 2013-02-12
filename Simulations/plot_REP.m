function [] = plot_REP(Big_phi, REP_cell,phi,MIP_cell, M, plot_panel)


% op_min = options(9);
op_min = 1;



N = size(phi,1);

% if op_figures
%     figure(20+fig_st);
%     plot_phi(phi,MIP_cell, M)
%     title('\phi  For Each Purview','FontSize',20);
%     xlabel('Purview','FontSize',16)
%     ylabel('\phi','FontSize',20)
% end


    

        

r = N;
c = 2;


set(plot_panel,'Position',[0.18076477404403243,0.0,0.7833140208574739,1.0]);
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

scaling = r/4;


if scaling > 1
      
    panel_pos = get(plot_panel,'Position');
    
    panel_height = panel_pos(4)*scaling;
    panel_pos(4) = panel_height;
    panel_pos(2) = 1 - panel_pos(4);
    
    set(plot_panel,'Position',panel_pos);
    
    %%% PICK UP HERE FIGURING OUT HOW TO SCALE THIS BITCH

    left_marg_col1 = .03;
    left_marg_col2 = .57;
    top_start = 1 - .45/panel_height;
    vert_marg = .1/panel_height;
    plot_h = .075/panel_height;
    plot_w = .4;
    
    setappdata(plot_panel,'limit',top_start-(vert_marg + plot_h)*(r-1)+.03)
    
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



%     Big_phi = sum(phi);
sy = ['Big Phi=',num2str(Big_phi)];
    
  assignin('base','MIP_Cell',MIP_cell);  
for i = 1:N
        

        
        prob = REP_cell{i,1};
        prob_prod = REP_cell{i,2};
        
        BR_w = prob{1};
        FR_w = prob{2};
        BR_p = prob_prod{1};
        FR_p = prob_prod{2};
        phi_b = phi(i,2);
        phi_f = phi(i,3);
        
        [string_p string] = make_title_fb(MIP_cell{i},0,1);
        
        s_title = cell(2,1);
        s_title_p = cell(2,1);
%         s_title{1}= [string{1},': \phi_b=',num2str(phi_b)];
%         s_title{2}= [string{2},': \phi_f=',num2str(phi_f)];
%         s_title_p{1} = [string_p{1},': \phi_b=',num2str(phi_b)];
%         s_title_p{2} = [string_p{2},': \phi_f=',num2str(phi_f)];
        
        s_title{1}= {['FULL:    ', string{1}, ',  \phi_b =',num2str(phi_b)],['PARTITIONED:    ', string_p{1}]};
        s_title{2}= {['FULL:    ', string{2}, ',  \phi_f =',num2str(phi_f)],['PARTITIONED:    ', string_p{2}]};
%          s_title{1} = {'test1','test2'};
%          s_title{2} = {'test3','test4'};
%         s = [string{3}, ', ', string_p{3},', phi = ',num2str(phi(i))];
%         fprintf('%s\n',s);
        

%             i_rev = N-i+1;
%             fig_pi = floor((i_rev-1)/fig_max);
%             fig_i = i_rev - fig_pi*fig_max;
            pos_BR = pos_vec(:,i,1)';
            pos_FR = pos_vec(:,i,2)';

%             figure(fig_st+fig_pi)
%             set(gcf,'Position',pos_fig)
            if mod(i,min(r,N)) == 0
                labelON = 1;
            else
                labelON = 0;
            end
            BR = ones(length(BR_w),2);
            BR(:,1) = BR_w; BR(:,2) = BR_p;
            FR = ones(length(FR_w),2);
            FR(:,1) = FR_w; FR(:,2) = FR_p;
            plot_BRFR(BR,FR,pos_BR,pos_FR,s_title,labelON,plot_panel)
    %         plot_BRFR(BR_w,FR_w,pos_BR',pos_FR',s_title,labelON)

            if i == N

                uicontrol('Style', 'text',...
                'String', {['Core Concept Distributions For ' mod_mat2str(M)],sy},... %replace something with the text you want
                'Units','normalized',...
                'FontSize',16,...
                'BackgroundColor','w',...
                'Position', [0.33 1-.2/panel_height 0.33 0.12/panel_height],...
                'Parent',plot_panel,...
                'Clipping','on');

    %             uicontrol('Style', 'text',...
    %             'String', 'PAST <-- CURRENT',... %replace something with the text you want
    %             'Units','normalized',...
    %             'FontSize',13,...
    %             'BackgroundColor','w',...
    %             'Position', [0.15 0.79 0.28 0.05]);

    %             uicontrol('Style', 'text',...
    %             'String', 'CURRENT --> FUTURE',... %replace something with the text you want
    %             'Units','normalized',...
    %             'FontSize',13,...
    %             'BackgroundColor','w',...
    %             'Position', [0.59 0.79 0.28 0.05]);

                my_legend = legend('Full Concept','Partitioned Concept');
                set(my_legend, 'Position',[.75 1-.2/panel_height .1 .18/panel_height],...
                'Parent',plot_panel,...
                'Clipping','on');

            end

            x0 = combine(MIP_cell{i}{1,2},MIP_cell{i}{2,2});
            subtit_pos = [.46 pos_BR(2) .06 .10/panel_height];

            uicontrol('Style', 'text',...
            'String', {mod_mat2str(x0),['phi = ', num2str(min(phi_f,phi_b))]},... %replace something with the text you want
            'Units','normalized',...
            'FontSize',13,...
            'BackgroundColor','w',...
            'Position', subtit_pos,...
            'Parent',plot_panel,...
            'Clipping','on');

end

end

                


function [] = plot_BRFR(BR,FR,pos_BR,pos_FR,s_title,labelON,plot_panel)


subplot('Position',pos_BR,'Parent',plot_panel)
h = bar(0:length(BR)-1,BR,'hist');
set(gca,'XTick',0:length(BR)-1)
set(gca,'YTickLabel',[0 .5 1])

states = convert(length(BR));
axis([-0.5 length(BR)-0.5 0 1.0])
title(s_title{1})
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


subplot('Position',pos_FR,'Parent',plot_panel);
h = bar(0:length(FR)-1,FR,'hist');
set(gca,'XTick',0:length(FR)-1)
set(gca,'YTickLabel',[0 .5 1])
states = convert(length(FR));
axis([-0.5 length(FR)-0.5 0 1.0])
title(s_title{2})
if labelON== 1
    set(gca,'XTickLabel',num2str(states,'%d')) % uncomment this to have a binary valued x-axis
    rotateXLabels( gca(), 90) % uncomment if binary values are used on the x-axis
%     xlabel('State (Node Order: [1...N])','FontSize',14)
%     xlab = xlabel('State (Node Order: [1...N])','FontSize',14,'Units','pixels');
%     set(xlab,'Position',get(xlab,'Position') - [0 100 0])

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
