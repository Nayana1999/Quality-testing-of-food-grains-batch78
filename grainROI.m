function [ IROI] = grainROI(I ,OI, cno)

ss=size(I);
sy=0;

left=0;
top=0;
right=0;
bottom=0;

%Locate Left Starting
for i=1:ss(1)
    temp(i)=0;
    for j=1:ss(2)
        if I(i,j)~=cno  
            temp(i)=temp(i)+1;
        else
            break;
        end
    end
end

left=min(temp);

%Locate top Starting
for j=1:ss(2)
    temp2(j)=0;
    for i=1:ss(1)
        if I(i,j)~=cno
            temp2(j)=temp2(j)+1;
        else
            break;
        end
        
    end
end

top=min(temp2);

%righSide End
for i=1:ss(1)
    temp(i)=0;
    for j=ss(2):-1:1
        if I(i,j)~=cno
            temp(i)=temp(i)+1;
        else
            break;
        end
    end
end
right=min(temp);


%Locate bottom Starting
for j=1:ss(2)
    temp2(j)=0;
    for i=ss(1):-1:1
        if I(i,j)~=cno 
            temp2(j)=temp2(j)+1;
        else
            break;
        end
        
    end
end

bottom=min(temp2);

IROI=imcrop(OI,[left+1 top+1  ss(2)-left-right ss(1)-top-bottom]);
%  figure
%  imshow(IROI);
%  title('Grain Segment');
end