function varargout = DehazeGUI(varargin)
% DEHAZEGUI MATLAB code for DehazeGUI.fig
%      DEHAZEGUI, by itself, creates a new DEHAZEGUI or raises the existing
%      singleton*.
%
%      H = DEHAZEGUI returns the handle to a new DEHAZEGUI or the handle to
%      the existing singleton*.
%
%      DEHAZEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEHAZEGUI.M with the given input arguments.
%
%      DEHAZEGUI('Property','Value',...) creates a new DEHAZEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DehazeGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DehazeGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DehazeGUI

% Last Modified by GUIDE v2.5 08-May-2016 11:48:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DehazeGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DehazeGUI_OutputFcn, ...
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


% --- Executes just before DehazeGUI is made visible.
function DehazeGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DehazeGUI (see VARARGIN)

% Choose default command line output for DehazeGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global imgInput;
global mask;

axes(handles.axes1);
imshow(rand(400,600));
axes(handles.axes2);
imshow(rand(200,300));


% UIWAIT makes DehazeGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DehazeGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global imgInput;
global mask;
[filename, pathname] = openImageFile();
if filename
    %imgInput = imread([pathname filename]);
    imgInput = imresize(imread([pathname filename]), 0.25);
    %imgInput = MaskingImage(imgInput);

    axes(handles.axes1);
    imshow(imgInput);
    title([filename]);
    
    siz = size(imgInput);
    siz = siz(1:2);
    center = floor(0.5*siz);

    % the mask is the region of pixels their distance the center is less ...
    mask = zeros(siz);
    mask(center(1), center(2)) = 1;
    mask = bwdist(mask, 'euclidean');
    mask = mask < 0.75*center(1) + 0.25*center(2);
    axes(handles.axes2);
    imshow(mask);
    
    imgHsv = rgb2hsv(imgInput);
    axes(handles.axes3);
    imshow(imgHsv(:,:,1));  title('Hue')
    axes(handles.axes4); 
    imshow(imgHsv(:,:,2));  title('Saturation');
    axes(handles.axes5);
    imshow(imgHsv(:,:,3));  title('Intensity Value');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global imgInput;

imgOutput = Dehaze(imgInput, 0.75);
axes(handles.axes2);
imshow(imgOutput/255);
title('Dehazed Image');

return;

if ~isempty(imgInput)
    [H, S, L] = rgb2hsl2(imgInput);
    s0 = 0.8;  % saturation level satruation_level = 0.8;
    h0 = 0.2; % the right bound of hue
    
    r = imgInput(:,:,1);
    g = imgInput(:,:,2);
    b = imgInput(:,:,3);
    
    Hue_rgb = zeros(3, 256);
    ROI = H < h0;
    ROI_Sat = ROI & S >= s0;
    ROI_Uns = ROI & S < s0;
    
    % the average r, g, b levels in saturated ROI
    % r_Sat_Mean = mean(mean(r(ROI_Sat)));
    %g_Sat_Mean = mean(mean(b(ROI_Sat)));
    %b_Sat_Mean = mean(mean(b(ROI_Sat)));
    
    %r_Uns_Mean = mean(mean(r(ROI_Uns)));
    %g_Uns_Mean = mean(mean(b(ROI_Uns)));
    %b_Uns_Mean = mean(mean(b(ROI_Uns)));
    
    % we divide hue interval [0, h0] into 256 bins
    % most of the code here can be vectorized
    H_bins = zeros(4, 256);
    
     [row, col] = size(H);
     for i = 1:row
        for j = 1:col
            if(ROI_Sat(i, j))
                binIndex = 1+floor(255*H(i, j));
                
                H_bins(1,binIndex) = H_bins(1, binIndex) + r(i, j);
                H_bins(2,binIndex) = H_bins(2, binIndex) + g(i, j);
                H_bins(3,binIndex) = H_bins(3, binIndex) + b(i, j);
                H_bins(4,binIndex) = H_bins(4, binIndex) + 1;     % pixel count
            end
        end
     end
     
    H_bins = H_bins + realmin('single');
    
    for i = 1:256
            H_bins(1, i) = H_bins(1, i)/H_bins(4, i);
            H_bins(2, i) = H_bins(2, i)/H_bins(4, i);
            H_bins(3, i) = H_bins(3, i)/H_bins(4, i);
    end
    R = r;
    G = g;
    B = b;
    for i = 1:row
        for j = 1:col
            if(ROI_Uns(i, j))
                lambda = S(i, j)/s0;
                binIndex = 1 + floor(255*H(i, j));
                R(i, j) = H_bins(1, binIndex) - lambda * ( H_bins(1, binIndex) - double(r(i, j)) );
                G(i, j) = H_bins(2, binIndex) - lambda * ( H_bins(2, binIndex) - double(g(i, j)) );
                B(i, j) = H_bins(3, binIndex) - lambda * ( H_bins(3, binIndex) - double(b(i, j)) );
                
%                 min_value = min(R(i, j),  G(i, j),  B(i, j));
%                 max_value = max(R(i, j),  G(i, j),  B(i, j));
            end
        end
    end
    
    imgOutput(:,:,1) = uint8(R);
    imgOutput(:,:,2) = uint8(G);
    imgOutput(:,:,3) = uint8(B);

end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global imgInput;

if ~isempty(imgInput)
    imgDehaze = Dehaze2(imgInput);
    axes(handles.axes2);
    imshow(imgDehaze);
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global imgInput;

if ~isempty(imgInput)
    axes(handles.axes2);
    title('Dehaze 3 rocessing ... ');
    pause(1);
    imgDehaze = Dehaze3(imgInput);
    imshow(imgDehaze);
    title('Dehazed');
    pause(1);
    

    imgEnhanced = AdjustLocalContrastTest(imgDehaze, 50);
    imshow(imgEnhanced);
    title('Dehaze + Local contrast adjustment');
    pause(1);
    
%     figure(99);
%     [h, s, v] = rgb2hsv(uint8(imgDehaze));
%     subplot(121); imshow(s); title('Saturation of dehazed image');
%     
%     figure(99);
%     [h, s, v] = rgb2hsv(uint8(imgEnhanced));
%     subplot(122); imshow(s); title('Saturation of enhanced image');
    
    return;
    
    for i = 1:3
        curChannel = imgDehaze(:,:,i);
        curChannel = cat(3, curChannel,curChannel,curChannel);
        axes(handles.axes4);
        imshow(imgDehaze);

        curChannel = AdjustLocalContrastTest(curChannel, 10);
        imgEnhanced(:,:,i) = curChannel(:,:,i);
    end
    imgEnhanced(:,:,1) = max(imgEnhanced(:,:,1), imgEnhanced(:,:,2));
    imgEnhanced(:,:,1) = max(imgEnhanced(:,:,1), imgEnhanced(:,:,3));
    imgEnhanced(:,:,3) = min(imgEnhanced(:,:,2), imgEnhanced(:,:,3));
    
    axes(handles.axes2);
    imshow(imgEnhanced);
    title('Haze reducing + Local Contrast Adjustment');
    figure(99); imshow(imgEnhanced);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global imgInput;
global mask;

if ~isempty(imgInput)
    imgDehaze = Dehaze5(imgInput, mask);
    axes(handles.axes2);
    imshow(imgDehaze/255);
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
