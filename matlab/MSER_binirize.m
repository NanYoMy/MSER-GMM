%mser binirize with differnet threshold
%??????
clear;
clc;
prepath='D:\mser\data';
image=load([prepath '\star.txt']);
maxT=max(max(image));
minT=min(min(image));
t=minT+(maxT-minT)*0.1
tmp=image;
imshow(uint8(image));

for i=1:size(tmp,1)
    for j=1:size(tmp,2)
        tmp(i,j)
        if tmp(i,j)>t
            tmp(i,j)=1;
        else
            tmp(i,j)=0;
        end
    end
end
imwrite(tmp,[prepath '\thread0.bmp']);
imshow(tmp);
 
tmp=image;
t=minT+(maxT-minT)*0.3
for i=1:size(tmp,1)
    for j=1:size(tmp,2)
        tmp(i,j)
        if tmp(i,j)>t
            tmp(i,j)=1;
        else
            tmp(i,j)=0;
        end
    end
end
imwrite(tmp,[prepath '\thread1.bmp']);
imshow(tmp); 

tmp=image;
t=minT+(maxT-minT)*0.5
for i=1:size(tmp,1)
    for j=1:size(tmp,2)
        tmp(i,j)
        if tmp(i,j)>t
            tmp(i,j)=1;
        else
            tmp(i,j)=0;
        end
    end
end
imwrite(tmp,[prepath '\thread2.bmp']);
imshow(tmp);
 
tmp=image;
t=minT+(maxT-minT)*0.9
for i=1:size(tmp,1)
    for j=1:size(tmp,2)
        tmp(i,j)
        if tmp(i,j)>t
            tmp(i,j)=1;
        else
            tmp(i,j)=0;
        end
    end
end
imwrite(tmp,[prepath '\thread3.bmp']);
imshow(tmp);

