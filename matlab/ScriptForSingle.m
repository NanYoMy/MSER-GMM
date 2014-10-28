function error= ScriptForSingle(X,Y,NM,NV)

prepath='D:\mser\data';
error=zeros(4,2);
x=X;
y=Y;
NoiseMean=NM;
NoiseVariance=NV;
% image=smallstary_one(prepath,'single',6.5,7.5,100,15);
image=smallstary_one(prepath,'single',x,y,NoiseMean,NoiseVariance);
biImage=binirization(prepath,'bi_single',image);
candidate=candidate_regions_extraction(prepath,'candidate',image,biImage);
system('D:\mser\code2\Release\code2.exe')
%run mser by code2.exe
%this method is used to call for the gmm2 to localize the position 

MULTIPLER=0.002;
samplesGenerator=load([prepath '\largetMSERRegion.txt']);
%samplesGenerator=binarize(samplesGenerator);
%imshow(samplesGenerator);
center=load([prepath '\highestIntensity.txt']);
%generater the sample
count=1;
sample=[];
for i=1:size(samplesGenerator)
    %这里开平方后就能得到与加权平方一样的效果
            for k=1:samplesGenerator(i,3)^2*MULTIPLER
                sample(count,:)=[samplesGenerator(i,2)-0.5 samplesGenerator(i,1)-0.5]; 
                count=count+1;
            end
%             for k=1:samplesGenerator(i,3)*MULTIPLER
%                 sample(count,:)=[samplesGenerator(i,2)-0.5 samplesGenerator(i,1)-0.5]; 
%                 count=count+1;
%             end
end
center=unique(center(:,1:2),'rows');
%center=2;
gmmMiu=gmm2(sample,center);
gmmX=gmmMiu(1,1);
gmmY=gmmMiu(1,2);
error(1,:)=[abs(gmmX-x) abs(gmmY-y)];

% fid=fopen([prepath '\gmmResult.txt'],'w+');
% for i=1:size(miu,1)
%     fprintf(fid,'%f %f \r\n',gmmMiu(i,1),gmmMiu(i,2));
% end
% fclose(fid);

%this method is use gaussian fitting
% clear;
% clc;
% prepath='D:\mser\data';
% samplesGenerator=load([prepath '\largetMSERRegion.txt']);
image=zeros(16,16);
for i=1:size(samplesGenerator)
    image(samplesGenerator(i,2),samplesGenerator(i,1))=samplesGenerator(i,3);
end
GaussianMu=GaussianFiting(image);
% fid=fopen([prepath '\guassianFit.txt'],'w+');
% fprintf(fid,'%f %f \r\n',mu(1),mu(2));
% fclose(fid);
%this method is use 
error(2,:)=[abs(GaussianMu(1)-x) abs(GaussianMu(2)-y)];
% clear;
% clc;
% prepath='D:\mser\data';
% samplesGenerator=load([prepath '\largetMSERRegion.txt']);
image=zeros(16,16);
for i=1:size(samplesGenerator)
    image(samplesGenerator(i,2),samplesGenerator(i,1))=samplesGenerator(i,3);
end

sum=0;
powerSum=0;
weightSumX=0;
weightSumY=0;
powerWSX=0;
powerWSY=0;

for i=1:16
    for j=1:16
        weightSumX=i*image(i,j)+weightSumX;
        powerWSX=i*image(i,j)^2+powerWSX;
        weightSumY=j*image(i,j)+weightSumY;
        powerWSY=j*image(i,j)^2+powerWSY;
        sum=image(i,j)+sum;
        powerSum=image(i,j)^2+powerSum;
        
    end
end

cx=weightSumX/sum-0.5;
cy=weightSumY/sum-0.5;
pcx=powerWSX/powerSum-0.5;
pcy=powerWSY/powerSum-0.5;
% fid=fopen([prepath '\cog.txt'],'w+');
% fprintf(fid,'%f %f \r\n',cx,cy);
% fclose(fid);
% fid=fopen([prepath '\powerCog.txt'],'w+');
% fprintf(fid,'%f %f \r\n',pcx,pcy);
% fclose(fid);
error(3,:)=[abs(cx-x) abs(cy-y)];
error(4,:)=[abs(pcx-x) abs(pcy-y)];
% clc;
% clear;
% prepath='D:\mser\data';
% %求解所有的
% a=load([prepath '\cog.txt']);
% b=load([prepath '\powerCog.txt']);
% c=load([prepath '\guassianFit.txt']);
% d=load([prepath '\gmmResult.txt']);

