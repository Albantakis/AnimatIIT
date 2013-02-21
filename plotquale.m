function varargout = plotquale(varargin)
% PLOTQUALE MATLAB code for plotquale.fig
%      PLOTQUALE, by itself, creates a new PLOTQUALE or raises the existing
%      singleton*.
%
%      H = PLOTQUALE returns the handle to a new PLOTQUALE or the handle to
%      the existing singleton*.
%
%      PLOTQUALE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTQUALE.M with the given input arguments.
%
%      PLOTQUALE('Property','Value',...) creates a new PLOTQUALE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotquale_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotquale_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotquale

% Last Modified by GUIDE v2.5 11-Jul-2012 17:42:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotquale_OpeningFcn, ...
                   'gui_OutputFcn',  @plotquale_OutputFcn, ...
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


% --- Executes just before plotquale is made visible.
function plotquale_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotquale (see VARARGIN)

% Choose default command line output for plotquale
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plotquale wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plotquale_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in concept_label_list.
function concept_label_list_Callback(hObject, eventdata, handles)
% hObject    handle to concept_label_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns concept_label_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from concept_label_list


% --- Executes during object creation, after setting all properties.
function concept_label_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to concept_label_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in upload_button.
function upload_button_Callback(hObject, eventdata, handles)
% hObject    handle to upload_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = uigetfile('*.mat');
load(filename,'concepts');
plotmatrix(concepts);
linkdata on
brush on
