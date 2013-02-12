%VERSION v0.04 - 11.14.2012

function varargout = iit(varargin)
% IIT MATLAB code for iit.fig
%      IIT, by itself, creates a new IIT or raises the existing
%      singleton*.
%
%      H = IIT returns the handle to a new IIT or the handle to
%      the existing singleton*.
%
%      IIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IIT.M with the given input arguments.
%
%      IIT('Property','Value',...) creates a new IIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iit

% Last Modified by GUIDE v2.5 11-Jan-2013 15:14:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iit_OpeningFcn, ...
                   'gui_OutputFcn',  @iit_OutputFcn, ...
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


% --- Executes just before iit is made visible.
function iit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iit (see VARARGIN)

% Choose default command line output for iit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set initial View
set(handles.Options,'Visible','Off')
set(handles.connectivity_mat,'Visible','Off')
set(handles.connect_mat_title,'Visible','Off')
set(handles.logic_types,'Visible','Off')
set(handles.logic_text,'Visible','Off')
set(handles.noise,'Visible','Off')
set(handles.noise_text,'Visible','Off')


% UIWAIT makes iit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function num_nodes_Callback(hObject, eventdata, handles)
% hObject    handle to num_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_nodes as text
%        str2double(get(hObject,'String')) returns contents of num_nodes as a double

if isnan(str2double(get(hObject,'String'))) || ~isposintscalar(str2double(get(hObject,'String')))

    set(handles.warning,'String','Number of nodes must be a positive integer.');
    set(hObject,'String',num2str(size(get(handles.TPM,'Data'),1)));
    
else
    set(handles.warning,'String','');
    tpm_choices = cellstr(get(handles.tpm_type_menu,'String'));
    tpm_choice = tpm_choices{get(handles.tpm_type_menu,'Value')};
    updateTPMview(handles, tpm_choice)
    updateCurrentStateView(handles)
    updateLogicTypesView(handles)
    updateConnectivityView(handles)
end



% --- Executes during object creation, after setting all properties.
function num_nodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in net_definition_method.
function net_definition_method_Callback(hObject, eventdata, handles)
% hObject    handle to net_definition_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns net_definition_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from net_definition_method

contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};

if strcmp(selection,'TPM / Connections')
    
    set(handles.TPM,'Enable','On')
    set(handles.logic_types,'Enable','Off')
    set(handles.noise,'Enable','Off')
    
elseif strcmp(selection,'Connections / Logic Mechanisms')
    
    user_response = confirm_def_switch('Title','Confirm Definition Change');
    
    switch lower(user_response)
        case 'no'
            return
        case 'yes'
            set(handles.TPM,'Enable','Off')
            set(handles.logic_types,'Enable','On')
            set(handles.noise,'Enable','On')
            update_tpm_from_connections_logic(handles)
            tpm_type_menu_Callback(handles.tpm_type_menu, eventdata, handles)
    end
    
    
end




% --- Executes during object creation, after setting all properties.
function net_definition_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to net_definition_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% redraw the TPM based on current settings
function updateTPMview(handles, tpm_choice)

nNodes = str2double(get(handles.num_nodes,'String'));

tpm_old = get(handles.TPM,'Data');
tpm_size_old = size(tpm_old,2);

if strcmp(tpm_choice,'State X State')
    
    % update TPM

    tpm_size_new = 2^nNodes;
    % increase size
    if tpm_size_new > tpm_size_old

        tpm_new = eye(tpm_size_new);
        tpm_new(1:size(tpm_old,1),1:size(tpm_old,2)) = tpm_old;

    % decrease size
    else

        % resize
        tpm_new = tpm_old(1:tpm_size_new,1:tpm_size_new);

    end

    set(handles.TPM,'Data',tpm_new);

    % rename cols and rows
    names = cell(1,tpm_size_new);
    for i = 1:tpm_size_new
        names{i} = dec2bin(i-1,nNodes);
    end
    set(handles.TPM,'ColumnName',names,'RowName',names);
    set(handles.TPM,'ColumnEditable',true(1,tpm_size_new));

elseif strcmp(tpm_choice,'State X Node')

    
    if nNodes > tpm_size_old
        
        tpm_new = zeros(2^nNodes,nNodes);
        tpm_new(1:2^tpm_size_old,1:tpm_size_old) = tpm_old;
    else
        tpm_new = tpm_old(1:2^nNodes,1:nNodes);
    end
    
    set(handles.TPM,'Data',tpm_new);
    
    row_names = cell(1,2^nNodes);
    for i = 1:2^nNodes
        row_names{i} = dec2bin(i-1,nNodes);
    end
    set(handles.TPM,'RowName',row_names);
    col_names = cell(1,nNodes);
    for i = 1:nNodes
        col_names{i} = num2str(i);
    end
    set(handles.TPM,'ColumnName',col_names);
    set(handles.TPM,'ColumnEditable',true(1,nNodes));
end



    
    

% redraw the current state based on new input
function updateCurrentStateView(handles)

nNodes = str2double(get(handles.num_nodes,'String'));

% update current state table

cur_state_old = get(handles.cur_state,'Data');
cur_state_size_old = length(cur_state_old);

% increase size
if nNodes > cur_state_size_old
    
    cur_state_new = zeros(1,nNodes);
    cur_state_new(1:cur_state_size_old) = cur_state_old;

% decrease size
else
    
    cur_state_new = cur_state_old(1:nNodes);
    
end

set(handles.cur_state,'Data',cur_state_new);
set(handles.cur_state,'ColumnEditable',true(1,nNodes));



% --- Executes on selection change in view_select.
function view_select_Callback(hObject, eventdata, handles)
% hObject    handle to view_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns view_select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from view_select

view_choices = cellstr(get(hObject,'String'));
% turn all views off
for i = 1:length(view_choices)
    
    this_view = view_choices{i}(view_choices{i} ~= ' ');
    eval(['set(handles.' this_view ',''Visible'',''Off'')'])
    
end

% turn selected view on
selection = view_choices{get(hObject,'Value')};
selection = selection(selection ~= ' ');
eval(['set(handles.' selection ',''Visible'',''On'')'])




% --- Executes during object creation, after setting all properties.
function view_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in parallel_option_menu.
function parallel_option_menu_Callback(hObject, eventdata, handles)
% hObject    handle to parallel_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parallel_option_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parallel_option_menu


% --- Executes during object creation, after setting all properties.
function parallel_option_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parallel_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in state_option_menu.
function state_option_menu_Callback(hObject, eventdata, handles)
% hObject    handle to state_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns state_option_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from state_option_menu

if get(hObject,'Value') == 1
    set(handles.cur_state,'Enable','on')
else
    set(handles.cur_state,'Enable','off')
end


% --- Executes during object creation, after setting all properties.
function state_option_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to state_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in complex_option_menu.
function complex_option_menu_Callback(hObject, eventdata, handles)
% hObject    handle to complex_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns complex_option_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from complex_option_menu


% --- Executes during object creation, after setting all properties.
function complex_option_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to complex_option_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% THIS WILL GET CHANGED WHEN OPTIONS ARE NAILED DOWN

op_big_phi = get(handles.big_phi_alg_menu,'Value') - 1;

op_normalize = get(handles.normalization_menu,'Value');
op_normalize_small_phi = 0;
op_normalize_big_phi = 0;
if op_normalize == 1 || op_normalize == 2
    op_normalize_small_phi = 1;
end
if op_normalize == 1 || op_normalize == 3
    op_normalize_big_phi = 1;
end

op_complex = get(handles.complex_option_menu,'Value');
if op_complex == 2
    op_complex = 0;
end

op_small_phi = get(handles.small_phi_func_menu,'Value') - 1;
op_ave = get(handles.state_option_menu,'Value') - 1;
op_parallel = get(handles.parallel_option_menu,'Value') - 1;
op_parfor = get(handles.parfor_option_menu,'Value');
op_strongconn = get(handles.StrongConn_option_menu,'Value') - 1;
op_removal = get(handles.RemoveNoise_option_menu,'Value') - 1;

%options = [3 1 2 1 1 0 0 1 1 0 op_big_phi 0 ...
%           op_normalize_big_phi op_normalize_small_phi op_complex op_small_phi op_ave op_parallel];
options = [op_parallel op_ave op_complex op_small_phi op_big_phi op_normalize_small_phi ...
           op_normalize_big_phi 0 op_parfor op_strongconn op_removal 1 1 0 0 1 1 0];
       

tpm_choices = cellstr(get(handles.tpm_type_menu,'String'));
tpm_choice = tpm_choices{get(handles.tpm_type_menu,'Value')};       

tpm = get(handles.TPM,'Data');
num_states = size(tpm,1);
num_nodes = str2double(get(handles.num_nodes,'String'));

if size(tpm,2) == num_states
    
    
    new_tpm = zeros(num_states,num_nodes);
    
    for i = 1:num_states
        for j = 1:num_nodes
            for k = 1:num_states
                
                state = dec2bin(k-1,num_nodes);
                
                if strcmp(state(j),'1')
                    new_tpm(i,j) = new_tpm(i,j) + tpm(i,k);
                end
            end
        end
    end
    
    tpm = new_tpm;
end

% permute the tpm so that the row order matches 
% THIS WILL GO AWAY SOON (hopefully)
permuted_tpm = zeros(size(tpm));

for i = 1:num_states
    
    permuted_row_index = trans10(flipud(trans2(i-1,num_nodes)));
    permuted_tpm(permuted_row_index,:) = tpm(i,:);
    
end

tpm = permuted_tpm;
    

current_state = get(handles.cur_state,'Data')';
noise = str2double(get(handles.noise,'String'));

connectivity_matrix = get(handles.connectivity_mat,'Data');

% setup node strucs and cpts for each node
%inputs = struct('num',{1 2},'name',{'A_p' 'B_p'},'num_states',{2 2},'state_names',{{'0' '1'}},'logic_type',{2 3})
if strcmp(tpm_choice, tpm_choices(3)) == 1
    logic_types = num2cell(get(handles.logic_types,'Data'));
else
    logic_types = cell(num_nodes,1);
end    
% init struct array
nodes(2*num_nodes) = struct('num',2*num_nodes,'name',[num2str(num_nodes) '_c'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types{num_nodes},'cpt',[],...
                            'num_sys_nodes',num_nodes,'input_nodes',[]);
                        
                        
% make past node structs                        
for i = 1:num_nodes
    
    nodes(i) = struct('num',i,'name',[num2str(i) '_p'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types{i},'cpt',[],...
                            'num_sys_nodes',num_nodes,'input_nodes',[]);
    
end

% make current node structs and their tpms
for i = 1:num_nodes
    
    nodes(num_nodes + i) = struct('num',num_nodes + i,'name',[num2str(i) '_c'],'num_states',2,...
                            'state_names',{{'0' '1'}},'logic_type',logic_types{i},'cpt',[],...
                            'num_sys_nodes',num_nodes,'input_nodes',[]);

	input_nodes = 1:num_nodes;
    input_nodes_indices = input_nodes(logical(connectivity_matrix(i,:)));
    nodes(num_nodes + i).input_nodes = input_nodes_indices;

    nodes(num_nodes + i).cpt = cpt_factory_tpm(nodes(num_nodes + i), input_nodes_indices, nodes, 2*num_nodes, tpm);
    
%     if any(nodes(num_nodes + i).cpt ~= test_cpt)
%         disp('error')
%     end
    
    
    
end

assignin('base','nodes',nodes)
explorer_handle = findall(0,'tag','iit_explorer');
delete(explorer_handle)

drawnow

iit_run(tpm,connectivity_matrix,current_state,noise,options,nodes);



% --- Executes on selection change in big_phi_alg_menu.
function big_phi_alg_menu_Callback(hObject, eventdata, handles)
% hObject    handle to big_phi_alg_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns big_phi_alg_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from big_phi_alg_menu


% --- Executes during object creation, after setting all properties.
function big_phi_alg_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to big_phi_alg_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in big_phi_func_menu.
function big_phi_func_menu_Callback(hObject, eventdata, handles)
% hObject    handle to big_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns big_phi_func_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from big_phi_func_menu


% --- Executes during object creation, after setting all properties.
function big_phi_func_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to big_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in small_phi_func_menu.
function small_phi_func_menu_Callback(hObject, eventdata, handles)
% hObject    handle to small_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns small_phi_func_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from small_phi_func_menu


% --- Executes during object creation, after setting all properties.
function small_phi_func_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to small_phi_func_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noise_Callback(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise as text
%        str2double(get(hObject,'String')) returns contents of noise as a double

if isnan(str2double(get(hObject,'String')))
    set(handles.warning,'String','Noise must be a real number in [0,.5].');
    set(hObject,'String','0');
else
    noise = str2double(get(hObject,'String'));
    if noise < 0
        set(handles.warning,'String','Noise must be in [0,.5].');
        set(hObject,'String','0');
    elseif noise > .5
        set(handles.warning,'String','Noise must be in [0,.5].');
        set(hObject,'String','.5');
    else
        set(handles.warning,'String','');
    end
end

update_tpm_from_connections_logic(handles)




% --- Executes during object creation, after setting all properties.
function noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in normalization_menu.
function normalization_menu_Callback(hObject, eventdata, handles)
% hObject    handle to normalization_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normalization_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normalization_menu


% --- Executes during object creation, after setting all properties.
function normalization_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normalization_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in TPM.
function TPM_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to TPM (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

tpm = get(hObject,'Data');
tpm_choices = cellstr(get(handles.tpm_type_menu,'String'));
tpm_choice = tpm_choices{get(handles.tpm_type_menu,'Value')};

if any(isnan(eventdata.NewData) || any(eventdata.NewData < 0) || any(eventdata.NewData > 1))
    
    set(handles.warning,'String','Entries in the TPM must be real numbers in [0,1].');
    tpm(eventdata.Indices(1),eventdata.Indices(2)) = eventdata.PreviousData;
    set(hObject,'Data',tpm);
    
elseif strcmp(tpm_choice,'State X State') && any(sum(tpm,2) ~= 1)
    
    set(handles.warning,'String','Rows in the TPM must sum to 1');
    
else
    
    set(handles.warning,'String','');
    
end

    


% --- Executes when entered data in editable cell(s) in cur_state.
function cur_state_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to cur_state (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

state_vec = get(hObject,'Data');
if (eventdata.NewData ~= 0 && eventdata.NewData ~= 1)
    
    set(handles.warning,'String','Node states can only be 0 or 1')
    state_vec(eventdata.Indices(2)) = eventdata.PreviousData;
    set(hObject,'Data',state_vec)
else
    
    set(handles.warning,'String','');
    
end


% --- Executes on button press in upload_tpm.
function upload_tpm_Callback(hObject, eventdata, handles)
% hObject    handle to upload_tpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = uigetfile('*.mat');

if (filename ~= 0)
    
    load(filename)

    if exist('tpm')
        if size(tpm,2) == size(tpm,1) && isposintscalar(log2(size(tpm,1)))

            num_nodes = log2(size(tpm,1));

            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.tpm_type_menu,'Value',1); % state x state
            set(handles.TPM,'Data',tpm);
            updateTPMview(handles,'State X State');
            

        elseif size(tpm,1) == 2^size(tpm,2)

            num_nodes = size(tpm,2);

            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.tpm_type_menu,'Value',2); % state x node
            set(handles.TPM,'Data',tpm);

            updateTPMview(handles,'State X Node');
        end

        updateCurrentStateView(handles)
        
        set(handles.warning,'String','');
    else
        set(handles.warning,'String','No variable named ''tpm'' in that data file.')
    end
    
    tpm_def = 1;
    set(handles.net_definition_method,'Value',tpm_def);
    net_definition_method_Callback(handles.net_definition_method, eventdata, handles);
    set(handles.num_nodes,'String',num_nodes)
    num_nodes_Callback(handles.num_nodes,eventdata,handles)
    
end
    

% --- Executes on selection change in tpm_type_menu.
function tpm_type_menu_Callback(hObject, eventdata, handles)
% hObject    handle to tpm_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tpm_type_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tpm_type_menu

tpm_choices = cellstr(get(hObject,'String'));
tpm_choice = tpm_choices{get(hObject,'Value')};

if strcmp(tpm_choice,'State X State') || strcmp(tpm_choice,'State X Node')
    
    set(handles.tpm_text,'Visible','on')
    set(handles.TPM,'Visible','on')
    set(handles.connect_mat_title,'Visible','off')
    set(handles.connectivity_mat,'Visible','off')
    set(handles.logic_types,'Visible','off')
    set(handles.logic_text,'Visible','off')
    set(handles.noise,'Visible','Off')
    set(handles.noise_text,'Visible','Off')
    
else
    
    set(handles.tpm_text,'Visible','off')
    set(handles.TPM,'Visible','off')
    set(handles.connect_mat_title,'Visible','on')
    set(handles.connectivity_mat,'Visible','on')
    set(handles.logic_types,'Visible','on')
    set(handles.logic_text,'Visible','on')   
    set(handles.noise,'Visible','on')
    set(handles.noise_text,'Visible','on')
end
    
    

tpm = get(handles.TPM,'Data');
num_nodes = str2double(get(handles.num_nodes,'String'));
num_states = 2^num_nodes;

% if we want state x state and we were in state x node
if strcmp(tpm_choice,'State X State') && size(tpm,2) == num_nodes
    
    new_tpm = ones(num_states);
    
    % create state x state tpm
    for i = 1:num_states
        for j = 1:num_states
            
            state = num2str(dec2bin(j-1,num_nodes));
            
            for k = 1:num_nodes
                if strcmp(state(k),'1')
                    new_tpm(i,j) = new_tpm(i,j) * tpm(i,k);
                else
                    new_tpm(i,j) = new_tpm(i,j) * (1 - tpm(i,k));
                end
            end
        end
    end
    
    set(handles.TPM,'Data',new_tpm)
    
% if we want state x node and we were in state x state    
elseif strcmp(tpm_choice,'State X Node') && size(tpm,2) == num_states
    
    new_tpm = zeros(num_states,num_nodes);
    
    for i = 1:num_states
        for j = 1:num_nodes
            for k = 1:num_states
                
                state = dec2bin(k-1,num_nodes);
                
                if strcmp(state(j),'1')
                    new_tpm(i,j) = new_tpm(i,j) + tpm(i,k);
                end
            end
        end
    end
    
    set(handles.TPM,'Data',new_tpm)
end

updateTPMview(handles, tpm_choice)



% --- Executes during object creation, after setting all properties.
function tpm_type_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tpm_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function NetworkDefinition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NetworkDefinition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes when entered data in editable cell(s) in connectivity_mat.
function connectivity_mat_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to connectivity_mat (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

connectivity_mat = get(hObject,'Data');
if (eventdata.NewData ~= 0 && eventdata.NewData ~= 1)
    
    set(handles.warning,'String','Entries in connectivity matrix can only be 0 or 1')
    connectivity_mat(eventdata.Indices(1),eventdata.Indices(2)) = eventdata.PreviousData;
    set(hObject,'Data',connectivity_mat)
else
    
    if strcmp(get(handles.logic_types,'Enable'),'on')
        update_tpm_from_connections_logic(handles)
    end

    set(handles.warning,'String','');
    
end

% update_tpm_from_connections_logic(handles)


% --- Executes on button press in upload_connectivity.
function upload_connectivity_Callback(hObject, eventdata, handles)
% hObject    handle to upload_connectivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = uigetfile('*.mat');

if size(get(handles.TPM,'Data'),2) == str2double(get(handles.num_nodes,'String'))
    
    tpm_choice = 'State X Node';
    
else
    
    tpm_choice = 'State X State';
end


if (filename ~= 0)
    
    load(filename)

    if exist('connectivity_mat')
        
        if size(connectivity_mat,1) == size(connectivity_mat,2)

            num_nodes = size(connectivity_mat,1);
            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.connectivity_mat,'Data',connectivity_mat)
            
            updateCurrentStateView(handles)
            updateTPMview(handles,tpm_choice)
            updateConnectivityView(handles);
            set(handles.warning,'String','');
            
%             tpm_def = 1;
%             set(handles.net_definition,'Value',tpm_def);
%             net_definition_method_Callback(hObject, eventdata, handles);
        
        else
            
            set(handles.warning,'Connectivity Matrix is not square')
        
        end
    
    else
        set(handles.warning,'String','No variable named ''connectivity_mat'' in that data file.')
    end
    
    logical_def = 1;
    set(handles.net_definition_method,'Value',logical_def);
    net_definition_method_Callback(hObject, eventdata, handles);
    
end

function updateConnectivityView(handles)

nNodes = str2double(get(handles.num_nodes,'String'));

connect_mat_old = get(handles.connectivity_mat,'Data');
connect_mat_size_old = size(connect_mat_old,1);

% increase size
if nNodes > connect_mat_size_old

    connect_mat_new = ones(nNodes);
    connect_mat_new(1:size(connect_mat_old,1),1:size(connect_mat_old,2)) = connect_mat_old;

% decrease size
else

    % resize
    connect_mat_new = connect_mat_old(1:nNodes,1:nNodes);

end

set(handles.connectivity_mat,'Data',connect_mat_new);
set(handles.connectivity_mat,'ColumnEditable',true(1,nNodes));

if strcmp(get(handles.logic_types,'Enable'),'on')
    update_tpm_from_connections_logic(handles)
end
    




% --- Executes when entered data in editable cell(s) in logic_types.
function logic_types_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to logic_types (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

logic_vec = get(hObject,'Data');
if (eventdata.NewData < 0 || eventdata.NewData > 30)
    
    set(handles.warning,'String','Logic Types must be between 1 and 29')
    logic_vec(eventdata.Indices(2)) = eventdata.PreviousData;
    set(hObject,'Data',logic_vec)
    
else
    
    set(handles.warning,'String','');
    update_tpm_from_connections_logic(handles)
    
end

function update_tpm_from_connections_logic(handles)

updateTPMview(handles, 'State X Node')
logic_types = get(handles.logic_types,'Data');
nNodes = str2double(get(handles.num_nodes,'String'));
new_tpm = zeros(2^nNodes,nNodes);
connection_mat = get(handles.connectivity_mat,'Data');
noise = str2double(get(handles.noise,'String'));

for k = 1:2^nNodes
    
    % we must flip the vector so that the states come in the desired order
    % however, we have to remember to permute when we send to the backend
    % of the code... see run function - soon we will standarize this across
    % both
    x0 = flipud(trans2(k-1,nNodes));
    for i = 1:nNodes
        i_vec = logical(connection_mat(i,:));
        input_vec = x0(i_vec);
        new_tpm(k,i) = logic_gates(input_vec,logic_types(i),noise);
    end
    
end

set(handles.TPM,'Data',new_tpm)

% this is a bit of hack since we set the new tpm as state X node :(
if(get(handles.tpm_type_menu,'Value') == 1)
     
    stateXnode_view = 2;
    set(handles.tpm_type_menu,'Value',stateXnode_view);
    
end

function updateLogicTypesView(handles)

nNodes = str2double(get(handles.num_nodes,'String'));

log_types_old = get(handles.logic_types,'Data');
log_types_size_old = length(log_types_old);

% increase size
if nNodes > log_types_size_old
    
    log_types_new = ones(1,nNodes);
    log_types_new(1:log_types_size_old) = log_types_old;

% decrease size
else
    
    log_types_new = log_types_old(1:nNodes);
    
end

set(handles.logic_types,'Data',log_types_new);
set(handles.logic_types,'ColumnEditable',true(1,nNodes));


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tpm = get(handles.TPM,'Data');
tpm_view = get(handles.tpm_type_menu,'Value');
logic_types = get(handles.logic_types,'Data');
noise = get(handles.noise,'String');
connectivity_mat = get(handles.connectivity_mat,'Data');
cur_state = get(handles.cur_state,'Data');
net_definition = get(handles.net_definition_method,'Value');

options_handles = findobj('Style','popupmenu','Parent',handles.Options);
num_options = length(options_handles);
options_tags = cell(num_options,1);
options_values = zeros(num_options,1);
for i = 1:length(options_handles)
    
    options_tags{i} = get(options_handles(i),'Tag');
    options_values(i) = get(options_handles(i),'Value');
    
end
    
    
variable_names = {'tpm','tpm_view','logic_types','noise','connectivity_mat','cur_state','options_tags','options_values','net_definition'};
uisave(variable_names);


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = uigetfile('*.mat');

if (filename ~= 0)
    
    load(filename)
    
    if exist('tpm')
        if size(tpm,2) == size(tpm,1) && isposintscalar(log2(size(tpm,1)))

            num_nodes = log2(size(tpm,1));

            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.tpm_type_menu,'Value',1); % state x state
            set(handles.TPM,'Data',tpm);
            updateTPMview(handles,'State X State');


        elseif size(tpm,1) == 2^size(tpm,2)

            num_nodes = size(tpm,2);

            set(handles.num_nodes,'String',num2str(num_nodes))
            set(handles.tpm_type_menu,'Value',2); % state x node
            set(handles.TPM,'Data',tpm);
            updateTPMview(handles,'State X Node');
        end
    end

        
    if exist('tpm_view')
        set(handles.tpm_type_menu,'Value',tpm_view);
    end
    
    if exist('logic_types')
        set(handles.logic_types,'Data',logic_types)
    end
    
    if exist('noise')
        set(handles.noise,'String',noise)
    end
    
    if exist('cur_state')
        set(handles.cur_state,'Data',cur_state);
    end
    
    if exist('connectivity_mat')

            set(handles.connectivity_mat,'Data',connectivity_mat)

    end 
    
    if exist('options_tags') && exist('options_values')
                
        options_handles = findobj('Style','popupmenu','Parent',handles.Options);
        
        for i = 1:length(options_handles)
            
            tag = get(options_handles(i),'Tag');
            [found, index] = ismember('options_tags', tag)
            
            if found
                set(options_handles(i),'Value',options_values(index))
            end

        end
    end
    
    if exist('net_definition')
        set(handles.net_definition_method,'Value',net_definition)
    end
    
    net_definition_method_Callback(hObject, eventdata, handles);
        
    updateCurrentStateView(handles)
    updateConnectivityView(handles);
    updateLogicTypesView(handles);
        
end


% --- Executes on key press with focus on big_phi_alg_menu and none of its controls.
function big_phi_alg_menu_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to big_phi_alg_menu (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
