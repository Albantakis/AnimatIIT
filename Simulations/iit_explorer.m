function varargout = iit_explorer(varargin)
% IIT_EXPLORER MATLAB code for iit_explorer.fig
%      IIT_EXPLORER, by itself, creates a new IIT_EXPLORER or raises the existing
%      singleton*.
%
%      H = IIT_EXPLORER returns the handle to a new IIT_EXPLORER or the handle to
%      the existing singleton*.
%
%      IIT_EXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IIT_EXPLORER.M with the given input arguments.
%
%      IIT_EXPLORER('Property','Value',...) creates a new IIT_EXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iit_explorer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iit_explorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iit_explorer

% Last Modified by GUIDE v2.5 31-Jan-2013 18:16:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iit_explorer_OpeningFcn, ...
                   'gui_OutputFcn',  @iit_explorer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before iit_explorer is made visible.
function iit_explorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iit_explorer (see VARARGIN)

% set(hObject, 'Renderer', 'painters')

% Choose default command line output for iit_explorer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iit_explorer wait for user response (see UIRESUME)
% uiwait(handles.iit_explorer);

% Set initial View
set(handles.Concepts,'Visible','Off')
set(handles.SystemPartitions,'Visible','Off')

% load data
if nargin == 4 && isstruct(varargin{1})
    handles.data = varargin{1};
    guidata(hObject,handles)
else
    fprintf('Please load a struct when opening IIT_EXPLORER\n');
    handles.output = 'ERROR';
    guidata(hObject,handles)
    close(gcf)
    return
end

handles.mip_axes = [];
guidata(hObject,handles);
handles.export_plot = 0;
guidata(hObject,handles);

num_states = 2^handles.data.num_nodes;
all_states = length(handles.data.Big_phi_M) > 1;

% setup state listbox
if all_states
    states = cell(1,num_states + 1);
    states{1} = 'Average';
    for i = 1:num_states
        states{i + 1} = dec2bin(i-1,handles.data.num_nodes);
    end
else
    states = {mod_mat2str(handles.data.current_state')};
    set(handles.state_list,'Enable','off')
end
set(handles.state_list,'String',states);

% setup subset listbox
nodes = cell(1,handles.data.num_nodes);
for i = 1:handles.data.num_nodes
    nodes{i} = num2str(i);
end
set(handles.nodes_list,'String',nodes)
set(handles.nodes_list,'Value',handles.data.Complex{1})

set(handles.overview_axes_panel,'Visible','off');
set(handles.summary_panel,'Visible','off');

refresh_subset_button_Callback(handles.refresh_subset_button,eventdata,handles)



% set(handles.overview_scroll_panel,'Parent',handles.overview_axes_panel)





% --- Outputs from this function are returned to the command line.
function varargout = iit_explorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = 'Explorer Closed.';


% --- Executes on selection change in view_menu.
function view_menu_Callback(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns view_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from view_menu

subset = handles.data.subset; subset_index = handles.data.subset_index;
state_index = handles.data.state_index;

view_choices = cellstr(get(hObject,'String'));
% turn all views off
for i = 1:length(view_choices)
    
    this_view = view_choices{i}(view_choices{i} ~= ' ');
    eval(['set(handles.' this_view ',''Visible'',''Off'')'])
    
end

% turn selected view on
view = view_choices{get(hObject,'Value')};
view_no_spaces = view(view ~= ' ');
eval(['set(handles.' view_no_spaces ',''Visible'',''On'')'])

if strcmp(view,'Overview')
    
    set(handles.overview_axes_panel,'Visible','off')
    
    drawnow
    
    current_display_elements = allchild(handles.overview_axes_panel);

    for i = 1:length(current_display_elements)
            delete(current_display_elements(i))
    end



    set(handles.overview_axes_text,'Visible','off')
    set(handles.big_phi_text,'String',['Big Phi = ' num2str(handles.data.Big_phi_M{state_index}(subset_index))])
    set(handles.big_phi_MIP_text,'String',['Big Phi MIP = ' num2str(handles.data.Big_phi_MIP{state_index}(subset_index))])
    set(handles.MIP_text,'String',{'MIP:',[mod_mat2str(handles.data.complex_MIP_M{state_index}{subset_index}) '-'...
                                            mod_mat2str(pick_rest(subset,handles.data.complex_MIP_M{state_index}{subset_index}))]})
    set(handles.sum_small_phi_text,'String',['Sum Small Phi = ' num2str(sum(handles.data.small_phi_M{state_index}{subset_index}(:,1)))])
    set(handles.num_core_concepts_text,'String',['# Core Concepts = ' num2str(sum(handles.data.small_phi_M{state_index}{subset_index}(:,1) ~= 0))])

    [IRR_REP IRR_phi IRR_MIP M_IRR] = IRR_points(handles.data.concepts_M{state_index},...
                                                 handles.data.small_phi_M{state_index},...
                                                 handles.data.concept_MIP_M{state_index},subset, subset_index);


                                           
    if handles.export_plot
        figure_handle = figure;
        panel = uipanel('Parent',figure_handle);
        set(handles.export_plot_button,'BackgroundColor',[0.9294    0.9294    0.9294]);
    else
        panel = handles.overview_axes_panel;
    end
    
    plot_REP(handles.data.Big_phi_M{state_index}(subset_index), IRR_REP, IRR_phi, IRR_MIP,...
                                        subset, panel)
	% reset export_panel flag
    handles.export_plot = 0;
    guidata(hObject,handles);


    set(handles.summary_panel,'Visible','on')
    set(handles.overview_axes_panel,'Visible','on')
    set(handles.panel_slider,'Value',1.0)
    
elseif strcmp(view,'Concepts')
    
    N = length(subset);
    
    M_cell = cell(2^N-1,1);
    
    k = 1;
    for i = 1:N 
        C = nchoosek(subset,i); 
        N_C = size(C,1);
        for j = 1:N_C % for all combos of size i
            x0 = C(j,:); % pick a combination
            M_cell{k} = x0;% store combo
            k = k + 1;
        end
    end

    num_concepts = length(M_cell);
    concept_names = cell(2*num_concepts,1);
    for i = 1:num_concepts
        concept_names{i} = [mod_mat2str(M_cell{i}) '_past'];
    end
    for i = num_concepts+1:2*num_concepts
        concept_names{i} = [mod_mat2str(M_cell{i-num_concepts}) '_future'];
    end
    set(handles.output_concepts_list,'String',concept_names)
    set(handles.output_concepts_list,'Value',[]);    

elseif strcmp(view,'System Partitions')
    
    partition_names = cell(length(handles.data.complex_MIP_all_M{state_index}{subset_index}),1);
    for i = 1:length(handles.data.complex_MIP_all_M{state_index}{subset_index})
        partition_names{i} = [mod_mat2str(handles.data.complex_MIP_all_M{state_index}{subset_index}{i}) '-'...
                             mod_mat2str(pick_rest(subset,handles.data.complex_MIP_all_M{state_index}{subset_index}{i}))];
    end
    set(handles.partition_list,'String',partition_names)
    
    select_MIP(handles, subset, state_index, subset_index)

    update_purviews_list(handles)
    
    plot_partition(handles);
    
end


% --- Executes during object creation, after setting all properties.
function view_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nodes_list.
function nodes_list_Callback(hObject, eventdata, handles)
% hObject    handle to nodes_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nodes_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nodes_list


% --- Executes during object creation, after setting all properties.
function nodes_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nodes_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in complex_button.
function complex_button_Callback(hObject, eventdata, handles)
% hObject    handle to complex_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)\

state_index = get_state_index(handles);
set(handles.nodes_list,'Value',handles.data.Complex{state_index})

refresh_subset_button_Callback(handles.refresh_subset_button, eventdata, handles)



function state_index = get_state_index(handles)

state_choice = get(handles.state_list,'Value');
if length(get(handles.state_list,'String')) == 1
    state_index = 1;
else
    state_index = trans10(flipud(trans2(state_choice - 2,handles.data.num_nodes)));
end


% % --------------------------------------------------------------------
% function brush_toggle_OffCallback(hObject, eventdata, handles)
% % hObject    handle to brush_toggle (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% brush off
% 
% 
% % --------------------------------------------------------------------
% function brush_toggle_OnCallback(hObject, eventdata, handles)
% % hObject    handle to brush_toggle (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% brush on
% 
% % below is an attempt to get 
% % for i = 1:length(handles.mip_axes)
% %     
% %     brushed = findall(handles.mip_axes{i},'tag','Brushing');
% %     set(handles.brushed,'Parent',handles.mip_plot_panel,'Clipping','on')
% %     
% % end


% --- Executes on selection change in partition_list.
function partition_list_Callback(hObject, eventdata, handles)
% hObject    handle to partition_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns partition_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from partition_list

update_purviews_list(handles)
set(handles.partition_plot_refresh,'BackgroundColor','g')


% --- Executes during object creation, after setting all properties.
function partition_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to partition_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in purviews_list.
function purviews_list_Callback(hObject, eventdata, handles)
% hObject    handle to purviews_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns purviews_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from purviews_list

set(handles.partition_plot_refresh,'BackgroundColor','g')


% --- Executes during object creation, after setting all properties.
function purviews_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to purviews_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_subset_button.
function refresh_subset_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_subset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.data.subset handles.data.subset_index handles.data.state_index] = get_subset_and_state(handles);
guidata(gcf,handles)

state_index = handles.data.state_index;

if isempty(handles.data.Complex{state_index})
    set(handles.overview_axes_panel,'Visible','off')
    set(handles.overview_axes_text,'String','This state in not realizable','Visible','on')
    return
end

view_menu_Callback(handles.view_menu, eventdata, handles)
    


    
function select_MIP(handles, subset, state_index, subset_index)

partition_names = get(handles.partition_list,'String');
MIP_string = [mod_mat2str(handles.data.complex_MIP_M{state_index}{subset_index}) '-'...
                                            mod_mat2str(pick_rest(subset,handles.data.complex_MIP_M{state_index}{subset_index}))];
MIP_index = find(strcmp(MIP_string,partition_names));
set(handles.partition_list,'Value',MIP_index)




%Larissa: This still needs to be checked if quali plot is correct!
function plot_partition(handles)

subset = handles.data.subset; subset_index = handles.data.subset_index;
state_index = handles.data.state_index;

[partition1 partition1_index partition2 partition2_index] = get_partitions(handles);

partition_index = 0;
for i = 1:length(partition1)-1
    partition_index = partition_index + nchoosek(length(subset),i);
end
[ignore_var additional] = intersect(nchoosek(subset,length(partition1)),partition1,'rows');
partition_index = partition_index + additional;

system_partition_text = ['System: ' mod_mat2str(subset)];
partition_choices = get(handles.partition_list,'String');
partition_text = partition_choices{get(handles.partition_list,'Value')};
system_partition_text = [system_partition_text ', Partition: ' partition_text];
big_phi_partition_text = ['Big Phi (partition) = ' num2str(handles.data.Big_phi_MIP_all_M{state_index}{subset_index}(partition_index,1))];

highlight_indices = get(handles.purviews_list,'Value');

if ~isempty(highlight_indices)
    purviews_list = get(handles.purviews_list,'String');
    selected_concepts_text = ['Selected Concepts: ' purviews_list{highlight_indices(1)}];
    for i = 2:length(highlight_indices)
        selected_concepts_text = [selected_concepts_text ', ' purviews_list{highlight_indices(i)}];
    end
else
    selected_concepts_text = 'Selected Concepts: None';
end

direction_choices = get(handles.past_future_list,'String');
direction_text = ['Temporal Direction: ' direction_choices{get(handles.past_future_list,'Value')}];

partition_info = {system_partition_text, big_phi_partition_text, selected_concepts_text, direction_text};

set(handles.partition_plot_info,'String',partition_info)

existing_axes = findobj(handles.mip_plot_panel,'Parent',handles.mip_plot_panel);
for k = 1:length(existing_axes)
   delete(existing_axes(k))
end

N = length(subset);

set(handles.mip_loading_text,'Visible','on')
set(handles.mip_plot_panel,'Visible','off')

drawnow

% get phi values for the whole
w_phi_all = handles.data.small_phi_M{state_index}{subset_index}(:,1)';
w_phi_concepts = w_phi_all(w_phi_all ~= 0);
%     IRR_whole = M_IRR_M{whole_i};


% get concepts for the whole
w_concept_dists_p = zeros(2^N,length(w_phi_concepts));
w_concept_dists_f = zeros(2^N,length(w_phi_concepts));

z = 1;
for i = 1:length(w_phi_all)
    if (w_phi_all(i) ~= 0)

        if ~isempty(handles.data.concepts_M{state_index}{subset_index,1}{i}{1})
            w_concept_dists_p(:,z) = handles.data.concepts_M{state_index}{subset_index,1}{i}{1};
        end
        if ~isempty(handles.data.concepts_M{state_index}{subset_index,1}{i}{2})
            w_concept_dists_f(:,z) = handles.data.concepts_M{state_index}{subset_index,1}{i}{2};
        end
        z = z + 1;
    end
end  

parts_phi_all = [handles.data.small_phi_M{state_index}{partition1_index}(:,1)' ...
                      handles.data.small_phi_M{state_index}{partition2_index}(:,1)'];

nIRR = sum(parts_phi_all ~= 0);

p_concept_dists_p = zeros(2^N,nIRR);
p_concept_dists_f = zeros(2^N,nIRR);
parts_phi_concepts = parts_phi_all(parts_phi_all ~= 0);

z = 1;
for k = 1:length(parts_phi_all)

    if (parts_phi_all(k) ~= 0)
        %we could change the if below to check against k instead...
        if(z <= sum(handles.data.small_phi_M{state_index}{partition1_index}(:,1) ~= 0))
            p_concept_dists_p(:,z) = expand_prob(handles.data.concepts_M{state_index}{partition1_index,1}{k}{1},subset,partition1);
            
            partition_rest = pick_rest(subset,partition1);
            fmaxent_rest = comp_pers_cpt(handles.data.network.nodes,[],partition_rest,[],'forward');                       
            p_concept_dists_f(:,z) = expand_prob_general(handles.data.concepts_M{state_index}{partition1_index,1}{k}{2},subset,partition1,fmaxent_rest(:));
        else
            k_offset = k - size(handles.data.small_phi_M{state_index}{partition1_index},1);
            p_concept_dists_p(:,z) = ...
                expand_prob(handles.data.concepts_M{state_index}{partition2_index,1}{k_offset}{1},subset,partition2);
            
            partition_rest = pick_rest(subset,partition2);
            fmaxent_rest = comp_pers_cpt(handles.data.network.nodes,[],partition_rest,[],'forward');                       
            p_concept_dists_f(:,z) = expand_prob_general(handles.data.concepts_M{state_index}{partition2_index,1}{k_offset}{2},subset,partition2,fmaxent_rest(:));
        end
        z = z + 1;

    end

end

plot_choices = get(handles.partition_plot_menu,'String');
plot_choice_index = get(handles.partition_plot_menu,'Value');
plot_choice = plot_choices{plot_choice_index};

if get(handles.past_future_list,'Value') == 1
    all_concepts = [w_concept_dists_p'; p_concept_dists_p'];
else
    all_concepts = [w_concept_dists_f'; p_concept_dists_f'];
end

part1_purviews = handles.data.purviews_M{state_index}{partition1_index};
part2_purviews = handles.data.purviews_M{state_index}{partition2_index};
n_part_purviews = length(part1_purviews)+length(part2_purviews);
part_purviews = cell(n_part_purviews,1);

for i = 1:n_part_purviews

    if i <= length(part1_purviews)
        part_purviews{i} = part1_purviews{i};
    else
        part_purviews{i} = part2_purviews{i - length(part1_purviews)};
    end

end

% options for plot view

% 3D & 2D Scatter - Variance
% 3D Scatter - Variance
% 3D Scatter - PCA (unavail)
% 2D Scatter - Variance
% Concept Bar Graphs

% reset panel position
set(handles.mip_plot_panel,'Position',[0.14600231749710313,0.01160541586073501,0.8122827346465816,0.9052224371373307]);

dim_choices = get(handles.state_selection_menu,'String');
dim_choice = dim_choices{get(handles.state_selection_menu,'Value')};

% display chosen plot view
if strcmp(plot_choice,'3D & 2D Scatter')

    if handles.export_plot
        figure_handle = figure;
        panel = uipanel('Parent',figure_handle);
        set(handles.export_plot_button,'BackgroundColor','white');
    else
        panel = handles.mip_plot_panel;
    end
    
    set(handles.partition_panel_slider,'Visible','off')
    conceptscatter3D2D(all_concepts,size(w_concept_dists_p,2), handles.data.purviews_M{state_index}{subset_index},...
            part_purviews, highlight_indices, panel, '2D3D', dim_choice);
        
	handles.export_plot = 0;
    guidata(handles.iit_explorer,handles);
    
elseif strcmp(plot_choice,'3D Scatter')
    
    
    if handles.export_plot
        figure_handle = figure;
        panel = uipanel('Parent',figure_handle);
        set(handles.export_plot_button,'BackgroundColor',[0.9294    0.9294    0.9294]);
    else
        panel = handles.mip_plot_panel;
    end

    
    set(handles.partition_panel_slider,'Visible','off')
    conceptscatter3D2D(all_concepts,size(w_concept_dists_p,2), handles.data.purviews_M{state_index}{subset_index},...
            part_purviews, highlight_indices, panel, '3D',dim_choice);
        
	handles.export_plot = 0;
    guidata(handles.iit_explorer,handles);
    
%     figure_handle = figure(999);
%     my_panel = uipanel('Parent',figure_handle);
%     set(handles.partition_panel_slider,'Visible','off')
%     conceptscatter3D2D(all_concepts,size(w_concept_dists_p,2), handles.data.purviews_M{state_index}{subset_index},...
%             part_purviews, highlight_indices, my_panel, '3D',dim_choice);    
        
        
elseif strcmp(plot_choice,'2D Scatter')

    
    if handles.export_plot
        figure_handle = figure;
        panel = uipanel('Parent',figure_handle);
        set(handles.export_plot_button,'BackgroundColor',[0.9294    0.9294    0.9294]);
    else
        panel = handles.mip_plot_panel;
    end
    
    set(handles.partition_panel_slider,'Visible','off')
    conceptscatter3D2D(all_concepts,size(w_concept_dists_p,2), handles.data.purviews_M{state_index}{subset_index},...
            part_purviews, highlight_indices, panel, '2D',dim_choice);  
        
	handles.export_plot = 0;
    guidata(handles.iit_explorer,handles);
    
elseif strcmp(plot_choice,'Concept Bar Graphs')
    
    if handles.export_plot
        figure_handle = figure;
        panel = uipanel('Parent',figure_handle);
        set(handles.export_plot_button,'BackgroundColor',[0.9294    0.9294    0.9294]);
    else
        panel = handles.mip_plot_panel;
    end    
    
    set(handles.partition_panel_slider,'Visible','on')
        
    plot_partition_bar(all_concepts,size(w_concept_dists_p,2), w_phi_concepts, parts_phi_concepts,...
        highlight_indices, panel, handles.data.purviews_M{state_index}{subset_index},...
            part_purviews, handles.data.Big_phi_MIP_all_M{state_index}{subset_index}(partition_index,1));
        
    set(handles.partition_panel_slider,'Value',1.0)
    
    handles.export_plot = 0;
    guidata(handles.iit_explorer,handles);

end


%     linkdata on

set(handles.mip_loading_text,'Visible','off')
set(handles.mip_plot_panel,'Visible','on') 

function [partition_p1 partition1_index partition_p2 partition2_index] = get_partitions(handles)

subset = handles.data.subset; subset_index = handles.data.subset_index;
state_index = handles.data.state_index;

partition_index = get(handles.partition_list,'Value');

partition_p1 = handles.data.complex_MIP_all_M{state_index}{subset_index}{partition_index};
partition1_index = convi(partition_p1) - 1;
partition_p2 = pick_rest(subset,partition_p1);
partition2_index = convi(partition_p2) - 1;


% --- Executes on selection change in state_list.
function state_list_Callback(hObject, eventdata, handles)
% hObject    handle to state_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns state_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from state_list


% --- Executes during object creation, after setting all properties.
function state_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to state_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function panel_slider_Callback(hObject, eventdata, handles)
% hObject    handle to panel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

shift = get(hObject,'Value');
position = get(handles.overview_axes_panel,'Position');
position(2) = (1 - position(4))*shift;
set(handles.overview_axes_panel,'Position',position)


% --- Executes during object creation, after setting all properties.
function panel_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in mip_button.
function mip_button_Callback(hObject, eventdata, handles)
% hObject    handle to mip_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subset = handles.data.subset; subset_index = handles.data.subset_index;
state_index = handles.data.state_index;

select_MIP(handles, subset, state_index, subset_index)

set(handles.partition_plot_refresh,'BackgroundColor','g')


% --- Executes on button press in partition_plot_refresh.
function partition_plot_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to partition_plot_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_partition(handles)
set(handles.partition_plot_refresh,'BackgroundColor',[.9294 .9294 .9294])


% --- Executes on button press in clear_purview_list.
function clear_purview_list_Callback(hObject, eventdata, handles)
% hObject    handle to clear_purview_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.purviews_list,'Value',[]);
% plot_partition(handles)

set(handles.partition_plot_refresh,'BackgroundColor','g')


% --- Executes on selection change in past_future_list.
function past_future_list_Callback(hObject, eventdata, handles)
% hObject    handle to past_future_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns past_future_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from past_future_list
set(handles.partition_plot_refresh,'BackgroundColor','g')

% --- Executes during object creation, after setting all properties.
function past_future_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to past_future_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in output_concepts_list.
function output_concepts_list_Callback(hObject, eventdata, handles)
% hObject    handle to output_concepts_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns output_concepts_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from output_concepts_list


% --- Executes during object creation, after setting all properties.
function output_concepts_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_concepts_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in output_concepts.
function output_concepts_Callback(hObject, eventdata, handles)
% hObject    handle to output_concepts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

concept_indices = get(handles.output_concepts_list,'Value');

subset = handles.data.subset; subset_index = handles.data.subset_index;
state_index = handles.data.state_index;

concept_count = length(get(handles.output_concepts_list,'String'))/2;

for i = 1:length(concept_indices)
    
    concept_names = cellstr(get(handles.output_concepts_list,'String'));
    name = concept_names{concept_indices(i)};
    name = ['purview_' mod_mat2str(subset) '_' name];
    name(name == '[') = [];
    name(name == ']') = [];
    name(name == ' ') = [];  
    if concept_indices(i) <= concept_count
        assignin('base',name,handles.data.concepts_M{state_index}{subset_index,1}{concept_indices(i)}{1})
    else
        assignin('base',name,handles.data.concepts_M{state_index}{subset_index,1}{concept_indices(i)-concept_count}{2})
    end

end


plot_partition(handles);

function [subset subset_index state_index] = get_subset_and_state(handles)

[subset subset_index] = get_subset(handles);

state_index = get_state_index(handles);

function [subset subset_index] = get_subset(handles)

subset = get(handles.nodes_list,'Value');
subset_index = convi(subset) - 1;

% make sure we don't return the empty set:
if subset_index == 0; subset_index = 1; end


% --- Executes on selection change in partition_plot_menu.
function partition_plot_menu_Callback(hObject, eventdata, handles)
% hObject    handle to partition_plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns partition_plot_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from partition_plot_menu

set(handles.partition_plot_refresh,'BackgroundColor','g')

if any(get(hObject,'Value') == [1 2 3])
    set(handles.state_selection_menu,'Visible','on')
else
    set(handles.state_selection_menu,'Visible','off')
end


% --- Executes during object creation, after setting all properties.
function partition_plot_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to partition_plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_purviews_list(handles)


subset_index = handles.data.subset_index;
state_index = handles.data.state_index;
[ignore_var partition1_index ignore_var partition2_index] = get_partitions(handles);

num_concepts_whole = length(handles.data.purviews_M{state_index}{subset_index,1});
num_concepts_part1 = length(handles.data.purviews_M{state_index}{partition1_index,1});
num_concepts_part2 = length(handles.data.purviews_M{state_index}{partition2_index,1});
num_concepts = num_concepts_whole + num_concepts_part1 + num_concepts_part2;
             
concept_names = cell(num_concepts,1);

for i = 1:num_concepts_whole
    concept_names{i} = [mod_mat2str(handles.data.purviews_M{state_index}{subset_index}{i}) '_whole'];
end

for i = num_concepts_whole + 1 : num_concepts_whole + num_concepts_part1
    concept_names{i} = [mod_mat2str(handles.data.purviews_M{state_index}{partition1_index}{i-num_concepts_whole}) '_part1'];
end

for i = num_concepts_whole + num_concepts_part1 + 1:num_concepts
    concept_names{i} = [mod_mat2str(handles.data.purviews_M{state_index}{partition2_index}{i-num_concepts_part1-num_concepts_whole}) '_part1'];
end

set(handles.purviews_list,'String',concept_names)
set(handles.purviews_list,'Value',[]);


% --- Executes on slider movement.
function partition_panel_slider_Callback(hObject, eventdata, handles)
% hObject    handle to partition_panel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

shift = get(hObject,'Value');
position = get(handles.mip_plot_panel,'Position');
position(2) = (1 - position(4))*shift;
set(handles.mip_plot_panel,'Position',position)


% --- Executes during object creation, after setting all properties.
function partition_panel_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to partition_panel_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in state_selection_menu.
function state_selection_menu_Callback(hObject, eventdata, handles)
% hObject    handle to state_selection_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns state_selection_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from state_selection_menu
set(handles.partition_plot_refresh,'BackgroundColor','g')


% --- Executes during object creation, after setting all properties.
function state_selection_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to state_selection_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in export_plot_button.
function export_plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to export_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.export_plot
    handles.export_plot = 0;
    set(hObject,'BackgroundColor',[.9294 .9294 .9294])
else
    handles.export_plot = 1;
    set(hObject,'BackgroundColor','red')
end
guidata(gcf,handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
