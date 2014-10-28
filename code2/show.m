%show the result
% pfx = fullfile('F:\Program Files\vlfeat-0.9.17','data','spots.jpg') ;
clc;
pfx = fullfile('D:\mser','data','StaryStaryNight256.bmp') ;
I = imread(pfx);
M=zeros(size(I));
result=load('D:\mser\code2\code2\posAndIntensity.txt');
[h,w]=size(result);
for i=1:1:h

    M(result(i,2),result(i,1))= M(result(i,2),result(i,1))+1; 
end


clf;imshow(I); hold on ; axis equal off; colormap gray ;
[c,h]=contour(M,(0:max(M(:)))+.5) ;
set(h,'color','y','linewidth',3) ;