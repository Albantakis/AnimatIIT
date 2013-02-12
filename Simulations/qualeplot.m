function varargout = qualeplot(concepts,labels,phis)
% tableplot.m -- Programmatic main function to set up a GUI
%                containing a uitable and an axes which
%                displays columns of the table as lines,
%                plus markers that show table selections.
%
% The following callbacks are also provided, as subfunctions:
%    plot1_callback   - Plot column selected on menu as line
%    select_callback  - Plot selected table data as markers 
%  Being subfunctions, they do not need handles passed to them.
%
%   Copyright 2008 The MathWorks, Inc.

% Create a figure that will have a uitable, axes and checkboxes
figure('Position', [100, 300, 600, 460],...
       'Name', 'TablePlot',...  % Title figure
       'NumberTitle', 'off',... % Do not show figure number
       'MenuBar', 'none');      % Hide standard menu bar menus

% Create a uitable on the left side of the figure
hlist = uicontrol('Style','listbox',...
                 'Units', 'normalized',...
                 'Position', [0.025 0.03 0.375 0.92],...
                 'String', labels,...
                 'ToolTipString',...
                 'Select concepts to highlight them on the plot',...
                 'CellSelectionCallback', {@select_callback});

% Create an axes on the right side; set x and y limits to the
% table value extremes, and format labels for the demo data.
[~,~,haxes] = plotmatrix(concepts);

title(haxes, 'Hourly Traffic Counts')   % Describe data set
% Prevent axes from clearing when new lines or markers are plotted
hold(haxes, 'all')

% Create an invisible marker plot of the data and save handles
% to the lineseries objects; use this to simulate data brushing.
hmkrs = plot(concepts, 'LineStyle', 'none',...
                    'Marker', 'o',...
                    'MarkerFaceColor', 'y',...
                    'HandleVisibility', 'off',...
                    'Visible', 'off');

               
% Subfuntions implementing the two callbacks
% ------------------------------------------

    function select_callback(hObject, eventdata)
    % hObject    Handle to uitable1 (see GCBO)
    % eventdata  Currently selected table indices
    % Callback to erase and replot markers, showing only those
    % corresponding to user-selected cells in table. 
    % Repeatedly called while user drags across cells of the uitable

        % hmkrs are handles to lines having markers only
        set(hmkrs, 'Visible', 'off') % turn them off to begin
        
        % Get the list of currently selected table cells
        sel = eventdata.Indices;     % Get selection indices (row, col)
                                     % Noncontiguous selections are ok
        selcols = unique(sel(:,2));  % Get all selected data col IDs
        table = get(hObject,'Data'); % Get copy of uitable data
        
        % Get vectors of x,y values for each column in the selection;
        for idx = 1:numel(selcols)
            col = selcols(idx);
            xvals = sel(:,1);
            xvals(sel(:,2) ~= col) = [];
            yvals = table(xvals, col)';
            % Create Z-vals = 1 in order to plot markers above lines
            zvals = col*ones(size(xvals));
            % Plot markers for xvals and yvals using a line object
            set(hmkrs(col), 'Visible', 'on',...
                            'XData', xvals,...
                            'YData', yvals,...
                            'ZData', zvals)
        end
    end
 end