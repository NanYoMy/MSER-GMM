%show the result
% pfx = fullfile('F:\Program Files\vlfeat-0.9.17','data','spots.jpg') ;
clc;
prepath='D:\mser\data\process2';
%pfx = fullfile('D:\mser','data','star.bmp') ;
I = imread([prepath '\double.bmp']);
M=zeros(size(I));
result=load([prepath '\posAndIntensity.txt']);
% result=load([prepath '\largetMSERRegion.txt']);
[h,w]=size(result);
for i=1:1:h

    M(result(i,2),result(i,1))= M(result(i,2),result(i,1))+1; 
end

clf;imshow(I); hold on ; axis equal off; colormap gray ;
[c,h]=contour(M,(0:max(M(:)))+.5) ;
set(h,'color','y','linewidth',3) ;