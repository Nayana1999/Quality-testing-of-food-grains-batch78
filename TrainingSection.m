%arranage training data set as matrix RxQ
clc
clear all
close all

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