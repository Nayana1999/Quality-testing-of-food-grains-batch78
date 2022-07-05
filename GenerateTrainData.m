function GenerateTrainData(TrainFolderPath)
J2=imread('bk2.jpg');

FolderLoaction = [TrainFolderPath 'BasmatiRice\'];
graintype = 1;%RICE
graingrade = 1;%BASMATI RICE
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
     % subplot(131),imshow(I);
      title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
      %subplot(132)
      %imshow(BW);
      title('graythresh')
      
      
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


     % subplot(133)
     % imshow(ColSeg)
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end


FolderLoaction = [TrainFolderPath 'RiceSonaMasuri\'];
graintype = 1;%RICE
graingrade = 2;%SonaMasuri RICE
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
      %subplot(131),imshow(I);
      %title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
     % subplot(132)
     % imshow(BW);
     % title('graythresh')
      
      
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


      %subplot(133)
     % imshow(ColSeg)
     % title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end



FolderLoaction = [TrainFolderPath 'GujratWheat\'];
graintype = 2;%Wheat
graingrade = 1;%GujratWheat
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
     % subplot(131),imshow(I);
     % title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
     % subplot(132)
    %  imshow(BW);
    %  title('graythresh')
      
      
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


     % subplot(133)
     % imshow(ColSeg)
    %  title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end



FolderLoaction = [TrainFolderPath 'KhapliWheat\'];
graintype = 2;%Wheat
graingrade = 2;%KhapliWheat
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
     % figure
     % subplot(131),imshow(I);
     % title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
     % subplot(132)
    %  imshow(BW);
    %  title('graythresh')
      
      
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


     % subplot(133)
     % imshow(ColSeg)
     % title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end



%Corn Data
FolderLoaction = [TrainFolderPath 'OrangeCorn\'];
graintype = 3;%Corn
graingrade = 1;%OrangeCorn
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
      %subplot(131),imshow(I);
     % title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
      %subplot(132)
      %%imshow(BW);
     % title('graythresh')
      
      
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


     % subplot(133)
     % imshow(ColSeg)
    %  title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end


FolderLoaction = [TrainFolderPath 'YellowCorn\'];
graintype = 3;%Corn
graingrade = 2;%YellowCorn
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
      %subplot(131),imshow(I);
      %title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
     % subplot(132)
     % imshow(BW);
     % title('graythresh')
      
      
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


    %  subplot(133)
     % imshow(ColSeg)
    %  title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end


%HorseGram

FolderLoaction = [TrainFolderPath 'HorseGramWhite\'];
graintype = 4;%HorseGram
graingrade = 2;%HorseGramWhite
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
      %subplot(131),imshow(I);
      %title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
    %  subplot(132)
    %  imshow(BW);
    %  title('graythresh')
      
      
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


     % subplot(133)
    %  imshow(ColSeg)
    %  title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end




FolderLoaction = [TrainFolderPath 'HorseGramBrown\'];
graintype = 4;%HorseGram
graingrade = 1;%HorseGramBrown
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
     % figure
     % subplot(131),imshow(I);
     % title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
    %  subplot(132)
    %  imshow(BW);
    %  title('graythresh')
      
      
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


     % subplot(133)
    %  imshow(ColSeg)
    %  title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end


FolderLoaction = [TrainFolderPath 'HorseGramWhite\'];
graintype = 4;%HorseGram
graingrade = 2;%HorseGramWhite
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
    %fullpath = strcat(pwd,FolderLoaction);
    lastletter = FileName(end);
    switch lastletter
        case 'A'
            grainqlty = 1;
        
        case 'B'
            grainqlty = 2;
       
        case 'C'
            grainqlty = 3;
        
        otherwise
            fprintf('Error, no such shape is found! Try again!\n')
    end
    image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc];
    SubFolderLocation = fullpath;
    disp(SubFolderLocation);
    SubFiles=dir(SubFolderLocation);
    for i=1:length(SubFiles)
      FileName=SubFiles(i).name;
      if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
      end
      fprintf('... %s %d%d%d\n', SubFiles(i).name,graintype,graingrade,grainqlty);
      J=imread([fullpath,SubFiles(i).name]);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
      %subplot(131),imshow(I);
      %title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
    %  subplot(132)
    %  imshow(BW);
     % title('graythresh')
      
      
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


      %subplot(133)
      %imshow(ColSeg)
     % title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);
    end
end


%IMPURITIES
FolderLoaction = [TrainFolderPath 'impurities\'];
graintype = 5;%Impurities
graingrade = 1;%Impurities
grainqlty = 1;%Impurities
Files=dir(FolderLoaction);
for k=1:length(Files)
    FileName=Files(k).name;
    
    if(strcmp(FileName,'.')||strcmp(FileName,'..')||strcmp(FileName,'Thumbs.db'))
       continue;
    end
    disp(FileName);
   image_loc = FileName
   % image_loc=strcat(FileName,'\');%BasmatiRiceA
    fullpath = [FolderLoaction,image_loc]
    
    J=imread(fullpath);
      OJ=J;
      J=J-J2;
      %figure,imshow(J);
     % title('GrainImage')
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
      %figure
     % subplot(131),imshow(I);
     % title('GrayImage')
      r=im2double(J(:,:,1));
        
      g=im2double(J(:,:,2));
      b=im2double(J(:,:,3));
      th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
      H=th;

      H(b>g)=2*pi-H(b>g);
      H=H/(2*pi);


      level=graythresh(H);

      BW=im2bw(I,level);
     % subplot(132)
     % imshow(BW);
    %  title('graythresh')
      
      
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


    %  subplot(133)
    %  imshow(ColSeg)
    %  title('segmented Color Image')
	  
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


		avgArea=0;
		avgmajaxis=0;
		avgminaxis=0;
		avgaspratio=0;
		avgrmean=0;
		avggmean=0;
		avgbmean=0; 

		%fprintf('GrainNo  Area  majoraxis minoraxis  aspratio  Rmean Gmean Bmean\n');
		fp=fopen('grainfinal.txt','a');
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

			avgArea=avgArea+Area;
			avgmajaxis=avgmajaxis+majaxis;
			avgminaxis=avgminaxis+minaxis;
			avgaspratio=avgaspratio+aspratio;
			avgrmean=avgrmean+rmean;
			avggmean=avggmean+gmean;
			avgbmean=avgbmean+bmean;

			Area=Area/10000;
			majaxis=majaxis/100;
			minaxis=minaxis/100;
			aspratio=aspratio/10;
			rmean=rmean/1000;
			gmean=gmean/1000;
			bmean=bmean/1000;

			
			fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

			fprintf(fp,'%d%d%d  ',graintype,graingrade,grainqlty);    
			fprintf(fp,'%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',Area,majaxis,minaxis,aspratio,rmean,gmean,bmean);    

		end

			avgArea=avgArea/graincount;
			avgmajaxis=avgmajaxis/graincount;
			avgminaxis=avgminaxis/graincount;
			avgaspratio=avgaspratio/graincount;
			avgrmean=avgrmean/graincount;
			avggmean=avggmean/graincount;
			avgbmean=avgbmean/graincount;

		%fprintf('Avg of features Extracted\n');
		%fprintf(1,'%d%d%d  ',graintype,graingrade,grainqlty);

		%fprintf('%8.4g %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',avgArea,avgmajaxis,avgminaxis,avgaspratio,avgrmean,avggmean,avgbmean);    

		fclose(fp);

end





end