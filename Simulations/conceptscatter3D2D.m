function [ax, height, extra_plots] = conceptscatter3D2D(x,nWholeConcepts, whole_purviews, part_purviews,...
                                            highlight_indices, parent_panel, view_option, dim_option)
% BASED ON GPLOTMATRIX



num_dims = min(size(x,2),8);
num_nodes = log2(size(x,2));

% princomp(x)

if strcmp(dim_option,'Variance')
    concept_value = var(x);
    
elseif strcmp(dim_option,'Mode')
%     concept_value = max(x);
    concept_value = sum(x,1);
end

[ignore_var state_ordering] = sort(concept_value,'descend');

nPartsConcepts = size(x,1) - nWholeConcepts;

whole = x(1:nWholeConcepts,:);
part = x(nWholeConcepts+1:end,:);

w_highlight_indices = highlight_indices(highlight_indices <= nWholeConcepts);
p_highlight_indices = highlight_indices(highlight_indices > nWholeConcepts) - nWholeConcepts;

w_nonhighlight_indices = setdiff(1:nWholeConcepts,w_highlight_indices);
p_nonhighlight_indices = setdiff(1:nPartsConcepts,p_highlight_indices);

whole_selected_labels = cell(length(w_highlight_indices),1);
for i = 1:length(w_highlight_indices)
    whole_selected_labels{i} = mod_mat2str(whole_purviews{w_highlight_indices(i)});
end

whole_nonselected_labels = cell(length(w_nonhighlight_indices),1);
for i = 1:length(w_nonhighlight_indices)
    whole_nonselected_labels{i} = mod_mat2str(whole_purviews{w_nonhighlight_indices(i)});
end

part_selected_labels = cell(length(p_highlight_indices),1);
for i = 1:length(p_highlight_indices)
    part_selected_labels{i} = mod_mat2str(part_purviews{p_highlight_indices(i)});
end

part_nonselected_labels = cell(length(p_highlight_indices),1);
for i = 1:length(p_nonhighlight_indices)
    part_nonselected_labels{i} = mod_mat2str(part_purviews{p_nonhighlight_indices(i)});
end

dims = size(x,2);

rows = 7; 
cols = rows;
extra_plots = rows - dims;

pos = get(parent_panel,'Position');
width = pos(3)/rows;
height = width;
space = .04; % 2 percent space between axes
xlim = repmat(cat(3,zeros(rows,1),ones(rows,1)),[rows 1 1]);
ylim = repmat(cat(3,zeros(rows,1)',ones(rows,1)'),[1 cols 1]);


x_bound = [0 1 0];
y_bound = [1 0 0];


ax = cell(nchoosek(num_dims,2)+1,1); % all pairs of dims plus the 3D plot
ax_index = 1;

if any(strcmp(view_option,{'2D','2D3D'}))
    
    for i = 8:-1:2 % count down from rows to 1
       for j = i-1:-1:1, % count down from cols to 1



            axPos = [(j-1)*width+space (rows-i+1)*height+space ...
                width*(1-space) height*(1-space)];
            ax{ax_index} = axes('Position',axPos, 'visible', 'on', 'Box','on','Parent',parent_panel,...
                'DrawMode','fast','Clipping','On');

            xlim(i,j,:) = get(ax{ax_index},'xlim');
            ylim(i,j,:) = get(ax{ax_index},'ylim');



            if (i <= num_dims)

                set(ax{ax_index},'Visible','on')
                state1 = state_ordering(i);
                state2 = state_ordering(j);


                plot(ax{ax_index},part(:,state2), ...
                    part(:,state1),'dr','Clipping','on');
%                 text(part(p_nonhighlight_indices,state2),part(p_nonhighlight_indices,state1),part_nonselected_labels)
                hold on;


                plot(ax{ax_index},whole(w_nonhighlight_indices,state2),...
                    whole(w_nonhighlight_indices,state1),'*b','Clipping','on')
                hold on;


                plot(ax{ax_index},whole(w_highlight_indices,state2), ...
                    whole(w_highlight_indices,state1),'og','MarkerSize',8,'Clipping','on');
                hold on;

                plot(ax{ax_index},part(p_highlight_indices,state2), ...
                    part(p_highlight_indices,state1),'og','MarkerSize',8,'Clipping','on');
                hold on;

                choices = nchoosek([1 2 3],2);

                for k = 1:size(choices,1)

                    hold on
                    plot(ax{ax_index},x_bound(choices(k,:)),y_bound(choices(k,:)),'k','Clipping','on');

                end
            else
                set(ax{ax_index},'Visible','off')
            end



            if j == 1 && i <= num_dims
                ylabel(ax{ax_index},dec2bin(state1-1,num_nodes))
            end
            if i == num_dims && j <= num_dims
                xlabel(ax{ax_index},dec2bin(state2-1,num_nodes))
            end

            set(ax{ax_index},'xlimmode','manual','ylimmode','manual','xgrid','off','ygrid','off',...
                    'xlim',[-.25 1.25],'ylim',[-.25 1.25],'xticklabel','','yticklabel','','Clipping','on')
            ax_index = ax_index + 1;

       end
    end
end


if any(strcmp(view_option,{'3D','2D3D'})) 
    
    if strcmp(view_option,'2D3D')
        row_pos = ceil(rows/2) - .5;
        col_pos = (rows - row_pos);
        size_scale = row_pos + 1.25;
        horiz_offset = 0;
    else
        row_pos = ceil(rows/8) - .5;
        col_pos = row_pos;
        size_scale = 24*row_pos;
        horiz_offset = .25;
    end

    axPos = [row_pos*width+space-horiz_offset col_pos*height+space ...
            width*(1-space)*size_scale height*(1-space)*size_scale];
    axes3D = axes('Position',axPos, 'visible', 'on', 'Box','on','Parent',parent_panel,'DrawMode','fast');

    ax{ax_index} = axes3D;


    scatter3(ax{ax_index},part(:,state_ordering(1)),part(:,state_ordering(2)),...
        part(:,state_ordering(3)),'Marker','d','MarkerEdgeColor','r','SizeData',75,'Clipping','on')
    hold on

    scatter3(ax{ax_index},whole(:,state_ordering(1)),whole(:,state_ordering(2)),...
        whole(:,state_ordering(3)),'Marker','*','MarkerEdgeColor','b','SizeData',75,'Clipping','on')

    hold on

    scatter3(ax{ax_index},whole(w_highlight_indices,state_ordering(1)),whole(w_highlight_indices,state_ordering(2)),...
        whole(w_highlight_indices,state_ordering(3)),'Marker','o','MarkerEdgeColor','g','SizeData',100,'Clipping','on')
    hold on

    scatter3(ax{ax_index},part(p_highlight_indices,state_ordering(1)),part(p_highlight_indices,state_ordering(2)),...
        part(p_highlight_indices,state_ordering(3)),'Marker','o','MarkerEdgeColor','g','SizeData',100,'Clipping','on')
    hold on

    xlabel(ax{ax_index},dec2bin(state_ordering(1)-1,num_nodes))
    ylabel(ax{ax_index},dec2bin(state_ordering(2)-1,num_nodes))
    zlabel(ax{ax_index},dec2bin(state_ordering(3)-1,num_nodes))


    set(ax{ax_index},'xlimmode','manual','ylimmode','manual',...
            'xlim',[-.25 1.25],'ylim',[-.25 1.25],'zlim',[-.25 1.25],...
            'CameraViewAngleMode','manual','Clipping','on')



    % plot tetrahedron bounds
    x_bound = [0 0 1 0];
    y_bound = [0 1 0 0];
    z_bound = [0 0 0 1];
    choices = nchoosek([1 2 3 4],2);

    for i = 1:size(choices,1)

        hold on
        plot3(ax{ax_index},x_bound(choices(i,:)),y_bound(choices(i,:)),z_bound(choices(i,:)),'k','Clipping','on');

    end
    if strcmp(view_option,'3D')
    %     text(part(p_nonhighlight_indices,state_ordering(1)),part(p_nonhighlight_indices,state_ordering(2)),...
    %         part(p_nonhighlight_indices,state_ordering(3)),part_nonselected_labels)
        text(part(p_highlight_indices,state_ordering(1))+.03,part(p_highlight_indices,state_ordering(2)),...
            part(p_highlight_indices,state_ordering(3)),part_selected_labels)
    %     text(whole(w_nonhighlight_indices,state_ordering(1)),whole(w_nonhighlight_indices,state_ordering(2)),...
    %         whole(w_nonhighlight_indices,state_ordering(3)),whole_nonselected_labels)
        text(whole(w_highlight_indices,state_ordering(1)),whole(w_highlight_indices,state_ordering(2)),...
            whole(w_highlight_indices,state_ordering(3)),whole_selected_labels)
    end
end
% linkdata on

% replace with real data
% whole = x(1:nWholeConcepts,:);
% part = x(nWholeConcepts+1:end,:);
% assignin('base','whole',whole)
% assignin('base','part',part)


% x(:) = in_data(:);
% assignin('base','x',x);


% ld = linkdata(gcf)
% fieldnames(ld)

% xlimmin = min(xlim(:,:,1),[],1); xlimmax = max(xlim(:,:,2),[],1);
% ylimmin = min(ylim(:,:,1),[],2); ylimmax = max(ylim(:,:,2),[],2);
% 
% % % Set all the limits of a row or column to be the same and leave 
% % % just a 5% gap between data and axes.
% inset = .05;
% for i=2:rows,
%   set(ax(i,1),'ylim',[ylimmin(i,1) ylimmax(i,1)])
%   dy = diff(get(ax(i,1),'ylim'))*inset;
%   set(ax(i,1:i-1),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
% end
% for j=1:cols-1,
%   set(ax(j+1,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
%   dx = diff(get(ax(1,j),'xlim'))*inset;
%   set(ax(j+1:rows,j),'xlim',[xlimmin(1,j)-dx xlimmax(1,j)+dx])
%   if ax2filled(j)
%      set(ax2(j),'xlim',[xlimmin(1,j)-dx xlimmax(1,j)+dx])
%   end
% end



% % Label plots one way or the other
% if (donames && ~isempty(xnam))
%    for j=1:cols
%       set(gcf,'CurrentAx',ax(j,j));
%       h = text((...
%           xlimmin(1,j)+xlimmax(1,j))/2, (ylimmin(j,1)+ylimmax(j,1))/2, -.1,...
%           xnam{j}, 'HorizontalAlignment','center',...
%           'VerticalAlignment','middle');
%    end
% else
%    if ~isempty(xnam)
%       for j=1:cols, xlabel(ax(rows,j),xnam{j}); end
%    end
%    if ~isempty(ynam)
%       for i=1:rows, ylabel(ax(i,1),ynam{i}); end
%    end
% end

% Ticks and labels on outer plots only
% set(ax(1:rows-1,:),'xticklabel','')
% set(ax(:,2:cols),'yticklabel','')
% set(BigAx,'XTick',get(ax(rows,1),'xtick'),'YTick',get(ax(rows,1),'ytick'), ...
%           'userdata',ax,'tag','PlotMatrixBigAx')

% Create legend if requested; base it on the top right plot
% if (doleg)
%    gn = gn(ismember(1:size(gn,1),g),:);
%    legend(ax(1,cols),gn);
% end

% Make BigAx the CurrentAxes
% set(gcf,'CurrentAx',BigAx)
% if ~hold_state,
%    set(gcf,'NextPlot','replace')
% end




% Also set Title and X/YLabel visibility to on and strings to empty
% set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
%  'String','','Visible','on')

% if nargout~=0,
%   h = hh;
%   if any(ax2filled)
%      ax = [ax; ax2(:)'];
%   end
% end

% -----------------------------
function datatipTxt = gplotmatrixDatatipCallback(obj,evt)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

dtcallbackdata = getappdata(target,'dtcallbackdata');
[BigAx,gnum,row,col] = dtcallbackdata{:};

ginds = getappdata(BigAx,'ginds');
xnam = getappdata(BigAx,'xnam');
ynam = getappdata(BigAx,'ynam');
xdat = getappdata(BigAx,'x');
ydat = getappdata(BigAx,'y');
XvsX = getappdata(BigAx,'XvsX');
gn = getappdata(BigAx,'gn');

gind = ginds{gnum};
obsind = gind(ind);

xvals = xdat(obsind,:);
yvals = ydat(obsind,:);

x = xvals(col);
y = yvals(row);

if x~=pos(1) || y~=pos(2)
    % Something is inconsistent, display default datatip.
    datatipTxt = {sprintf('X: %s',num2str(pos(1))),sprintf('Y: %s',num2str(pos(2)))};
else
    if isempty(xnam)
        xnam = cell(size(xdat,2),1);
        for i = 1:size(xdat,2)
            xnam{i} = sprintf('xvar%s',num2str(i));
        end
    end
    if isempty(ynam)
        ynam = cell(size(ydat,2),1);
        for i = 1:size(ydat,2)
            ynam{i} = sprintf('yvar%s',num2str(i));
        end
    end

    % Generate datatip text.
    datatipTxt = {
        [xnam{col},': ',num2str(x)],...
        [ynam{row},': ',num2str(y)],...
        '',...
        sprintf('Observation: %s',num2str(obsind)),...
        };

    if ~isempty(gn)
        datatipTxt{end+1} = ['Group: ',gn{gnum}];
    end
    datatipTxt{end+1} = '';

    xnamTxt = cell(length(xvals),1);
    for i=1:length(xvals)
        xnamTxt{i} = [xnam{i} ': ' num2str(xvals(i))];
    end
    datatipTxt = {datatipTxt{:}, xnamTxt{:}};
    
    if ~XvsX
        ynamTxt = cell(length(yvals),1);
        for i=1:length(yvals)
            ynamTxt{i} = [ynam{i} ': ' num2str(yvals(i))];
        end
        datatipTxt = {datatipTxt{:}, ynamTxt{:}};
    end

end

function [ogroup,glabel,gname,multigroup] = mgrp2idx(group,rows,sep); 
%MGRP2IDX Convert multiple grouping variables to index vector 
%   [OGROUP,GLABEL,GNAME,MULTIGROUP] = MGRP2IDX(GROUP,ROWS) takes 
%   the inputs GROUP, ROWS, and SEP.  GROUP is a grouping variable (numeric 
%   vector, string matrix, or cell array of strings) or a cell array 
%   of grouping variables.  ROWS is the number of observations. 
%   SEP is a separator for the grouping variable values. 
% 
%   The output OGROUP is a vector of group indices.  GLABEL is a cell 
%   array of group labels, each label consisting of the values of the 
%   various grouping variables separated by the characters in SEP. 
%   GNAME is a cell array containing one column per grouping variable 
%   and one row for each distinct combination of grouping variable 
%   values.  MULTIGROUP is 1 if there are multiple grouping variables 
%   or 0 if there are not. 
 
%   Tom Lane, 12-17-99 
%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 1.4 $  $Date: 2002/02/04 19:25:44 $ 
 
multigroup = (iscell(group) & size(group,1)==1); 
if (~multigroup) 
   [ogroup,gname] = grp2idx(group); 
   glabel = gname; 
else 
   % Group according to each distinct combination of grouping variables 
   ngrps = size(group,2); 
   grpmat = zeros(rows,ngrps); 
   namemat = cell(1,ngrps); 
    
   % Get integer codes and names for each grouping variable 
   for j=1:ngrps 
      [g,gn] = grp2idx(group{1,j}); 
      grpmat(:,j) = g; 
      namemat{1,j} = gn; 
   end 
    
   % Find all unique combinations 
   [urows,ui,uj] = unique(grpmat,'rows'); 
    
   % Create a cell array, one col for each grouping variable value 
   % and one row for each observation 
   ogroup = uj; 
   gname = cell(size(urows)); 
   for j=1:ngrps 
      gn = namemat{1,j}; 
      gname(:,j) = gn(urows(:,j)); 
   end 
    
   % Create another cell array of multi-line texts to use as labels 
   glabel = cell(size(gname,1),1); 
   if (nargin > 2) 
      nl = sprintf(sep); 
   else 
      nl = sprintf('\n'); 
   end 
   fmt = sprintf('%%s%s',nl); 
   lnl = length(fmt)-3;        % one less than the length of nl 
   for j=1:length(glabel) 
      gn = sprintf(fmt, gname{j,:}); 
      gn(end-lnl:end) = []; 
      glabel{j,1} = gn; 
   end 
end 

