clc;
close all;
clear all;

load Network;

[filename pathname]=uigetfile('.jpg','Select the Grain Image File For Testing');

fullpath=strcat(pathname,filename);

J=imread(fullpath);
J2=imread('bk2.jpg');
OJ=J;
J=J-J2;

figure
imshow(J);
title('Input GrainImage')

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

figure
imshow(BW)
title('segmented BW Image')

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
        
    
figure
%imshow(ColSeg)
title('segmented Color Image')

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

for i=1:2
    if pg(i)==max(pg)
        gq=i;
        break;
    end
end



disp('For Quality Prediction');
pg

qltfact=pg(gq)/sum(pg)*100;

m1='';
for i=1:length(gtype)
    if gtype(i)==5
        %msgbox('With IMPURITY');  
        m1='  With IMPURITY';
        break;
    end
end


m2=strcat(gt,' Contains ',num2str(qltfact),'%   A Quality grains  ',m1);

msgbox(m2);
% if qlty==1
%     msgbox('A-Grade');  
% else
%     msgbox('B-Grade');  
% end

 
 

pg

for i=1:5
    if pg(i)==max(pg)
        qlty=i;
        break;
    end
end