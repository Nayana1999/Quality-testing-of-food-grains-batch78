function varargout = FOODGRAINGUI(varargin)
% FOODGRAINGUI MATLAB code for FOODGRAINGUI.fig
%      FOODGRAINGUI, by itself, creates a new FOODGRAINGUI or raises the existing
%      singleton*.
%
%      H = FOODGRAINGUI returns the handle to a new FOODGRAINGUI or the handle to
%      the existing singleton*.
%
%      FOODGRAINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOODGRAINGUI.M with the given input arguments.
%
%      FOODGRAINGUI('Property','Value',...) creates a new FOODGRAINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FOODGRAINGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FOODGRAINGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FOODGRAINGUI

% Last Modified by GUIDE v2.5 20-Jun-2022 19:58:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FOODGRAINGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FOODGRAINGUI_OutputFcn, ...
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


% --- Executes just before FOODGRAINGUI is made visible.
function FOODGRAINGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FOODGRAINGUI (see VARARGIN)

% Choose default command line output for FOODGRAINGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FOODGRAINGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FOODGRAINGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_train_path.
function btn_train_path_Callback(hObject, eventdata, handles)
% hObject    handle to btn_train_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename path]= uigetfile('*.png;*.jpeg;*.jpg','Select Training Folder');

TrainFolderPath = path;

disp(TrainFolderPath);

save processed_data TrainFolderPath

% --- Executes on button press in btn_generate_train_data.
function btn_generate_train_data_Callback(hObject, eventdata, handles)
% hObject    handle to btn_generate_train_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('processed_data.mat')
fclose('all')
delete('grainfinal.txt');
GenerateTrainData(TrainFolderPath);
msgbox('Training Data Generated Successfully');

% --- Executes on button press in btn_perform_training.
function btn_perform_training_Callback(hObject, eventdata, handles)
% hObject    handle to btn_perform_training (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fp=fopen('grainfinal.txt','r');
Ptemp = fscanf(fp, '%g %g %g %g %g %g %g %g', [8 inf])
Ptemp=Ptemp';

msize=size(Ptemp);

for i=1:msize(1)
   P(i,1)=Ptemp(i,2);
   P(i,2)=Ptemp(i,3);
   P(i,3)=Ptemp(i,4);
   P(i,4)=Ptemp(i,5);
   P(i,5)=Ptemp(i,6);
   P(i,6)=Ptemp(i,7);
   P(i,7)=Ptemp(i,8);
end    
 
P=P';

for i=1:msize(1)
    Tc(i)=Ptemp(i,1);
end

T = ind2vec(Tc);

palmnet = newpnn(P,T,0.01);

save Network palmnet;

fprintf(1,'Training Over\n\n');
fclose(fp);
msgbox('Training Successfully Done');

% --- Executes on button press in btn_test_image.
function btn_test_image_Callback(hObject, eventdata, handles)
% hObject    handle to btn_test_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[filename pathname]=uigetfile('.jpg','Select the Grain Image File For Testing');

testfullpath=strcat(pathname,filename);

J=imread(testfullpath);
axes(handles.axes1)
imshow(J);

save processed_data testfullpath



% --- Executes on button press in btn_classify.
function btn_classify_Callback(hObject, eventdata, handles)
% hObject    handle to btn_classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load Network;
load processed_data

J=imread(testfullpath);
axes(handles.axes1)
imshow(J);

J2=imread('bk2.jpg');
OJ=J;
J=J-J2;

axes(handles.axes2)
imshow(J);
%title('Input GrainImage')

gsz=size(J);

I=rgb2gray(J);

for i=1:gsz(1)
    for j=1:gsz(2)
        if  I(i,j)<30
            I(i,j)=0;
        else
            I(i,j)=255;
        end
        
    end
end

r=im2double(J(:,:,1));
g=im2double(J(:,:,2));
b=im2double(J(:,:,3));

th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
H=th;

H(b>g)=2*pi-H(b>g);
H=H/(2*pi);


level=graythresh(H);


BW=im2bw(I,level);

axes(handles.axes3)
imshow(BW)

ColSeg=OJ;
for i=1:size(BW,1)
    for j=1:size(BW,2)
        if BW(i,j)==0
            ColSeg(i,j,1)=0;
            ColSeg(i,j,2)=0;
            ColSeg(i,j,3)=0;
        end
    
    end
end
        
    
axes(handles.axes4)
imshow(ColSeg)


[L num]=bwlabel(BW);

gainpatch=zeros(1,max(max(L)));

for i=1:gsz(1)    
    for j=1:gsz(2)    
    if L(i,j)>0
        gainpatch(L(i,j))=gainpatch(L(i,j))+1;
    end
    end
end

graincount=1;
for i=1:max(max(L))
    if gainpatch(i)>=70
        grainsegment(graincount)=i;
        graincount=graincount+1;
    end
end

graincount=graincount-1;
fprintf('Number of grains Identified %d\n',graincount);

for i=1:graincount
    ROI=grainROI(L,ColSeg , grainsegment(i));
    
    Area=gainpatch(grainsegment(i));
    
     if size(ROI,1)>size(ROI,2)
           majaxis=size(ROI,1);
           minaxis=size(ROI,2);
    else
           majaxis=size(ROI,2);
           minaxis=size(ROI,1);
     end
    
    aspratio=majaxis/minaxis;
    
    rmean=mean(mean(ROI(:,:,1)));
    gmean=mean(mean(ROI(:,:,2)));
    bmean=mean(mean(ROI(:,:,3)));

    Area=Area/10000;
    majaxis=majaxis/100;
    minaxis=minaxis/100;
    aspratio=aspratio/10;
    rmean=rmean/1000;
    gmean=gmean/1000;
    bmean=bmean/1000;

    
    fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    


te=[Area majaxis minaxis aspratio rmean gmean bmean];

Ptest=te;
Ptest=Ptest';

Y = sim(palmnet,Ptest);%Y = 121
pn2 = vec2ind(Y);

gtype(i)=fix(pn2/100);

ggrd(i)= fix(mod(pn2,100)/10);

gqlty(i)= fix(mod(mod(pn2,100),10));

end

%For finding the which Grain
pg=zeros(1,5);
for i=1:length(gtype)
    pg(gtype(i))=pg(gtype(i))+1;
end

pg

for i=1:5
    if pg(i)==max(pg)
        pn=i;%pn=1
        break;
    end
end

pn

%For Finding the grain grade
pg=zeros(1,3);
for i=1:length(ggrd)
    pg(ggrd(i))=pg(ggrd(i))+1;
end

pg

for i=1:2
    if pg(i)==max(pg)
        gg=i;
        break;
    end
end

gg


gt='';

 if pn==1%Rice
        if gg==1
            %msgbox('Basamati  Rice'); 
            gt='Basamati  Rice'; 
        else
            %msgbox('Sona Masuri  Rice');  
            gt='Sona Masuri Rice';  
        end
        
 elseif pn==2
   
     if gg==1
       %msgbox('Orange Corn');
       gt='Gujrat Wheat';
     else
       %msgbox('Yellow Corn'); 
       gt='Kapli Wheat';
     end
     
 elseif pn==3
         
     if gg==1
         %msgbox('GujratWheat'); 
         gt='Orange Corn';
     else
         %msgbox('KhapliWheat');
         gt='Yellow Corn';
     end
 
 elseif pn==4
     if gg==1
         gt='HorseGramBrown';
       %msgbox('HorseGramBrown'); 
     else
         gt='HorseGramWhite';
       %msgbox('HorseGramWhite');  
     end
 else
     msgbox('No match');
 end
 
%FOR QUALITY PREDICTION
pg=zeros(1,3);
for i=1:length(gqlty)
    pg(gqlty(i))=pg(gqlty(i))+1;
end

for i=1:3
    if pg(i)==max(pg)
        qlty=i;
        break;
    end
end

disp('For Quality Prediction');
qlty

qltfact=pg(1)/sum(pg)*100;

m1='';
for i=1:length(gtype)
    if gtype(i)==5
        %msgbox('With IMPURITY');  
        m1='  With IMPURITY';
        break;
    end
end

if qlty==1
    msgbox('A-Grade'); 
    gqty = 'A-Grade';	
else
    msgbox('B-Grade'); 
    gqty = 'B-Grade';	
end



m2=strcat(gt,' Contains  ',num2str(qltfact),'%  ',gqty,' Quality grains  ',m1);

set(handles.text10,'String',m2);
%msgbox(m2);


 
 

pg

for i=1:5
    if pg(i)==max(pg)
        qlty=i;
        break;
    end
end

% --- Executes on button press in btn_exit.
function btn_exit_Callback(hObject, eventdata, handles)
% hObject    handle to btn_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
