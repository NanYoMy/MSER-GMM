function error= ScriptForDouble(X1,Y1,X2,Y2,NM,NV)

prepath='D:\mser\data';

x1=X1;
y1=Y1;
x2=X2;
y2=Y2;
NoiseMean=NM;
NoiseVariance=NV;
% image=smallstary_one(prepath,'single',6.5,7.5,100,15);
image=smallstary_two(prepath,'double',x1,y1,x2,y2,NoiseMean,NoiseVariance);
biImage=binirization(prepath,'bi_double',image);
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
error=zeros(2,2);
error(1,1)=abs(gmmMiu(1,1)-x1);
error(1,2)=abs(gmmMiu(1,2)-y1);
error(2,1)=abs(gmmMiu(2,1)-x2);
error(2,2)=abs(gmmMiu(2,2)-y2);
