%高斯拟合，假设X与Y的方差一样
% clear;
% clc;
% prepath='D:\mser\data';
% image=load([prepath '\star.txt']);

function mu=GaussianFiting(image)
[w,h]=size(image);
A=zeros(h*w,4);
Y=zeros(h*w,1);
count=1;
for i=1:w
    for j=1:h
        
        if image(i,j)>0
            A(count,:)=[i^2+j^2 ,i , j , 1 ];
            tmp=image(i,j);
            Y(count)=log(tmp);
            count=count+1;
        end
    end
end

x=inv(A'*A)*A'*Y;
cx=-x(2)/2/x(1);
cx=cx-0.5;
cy=-x(3)/2/x(1);
cy=cy-0.5;
mu(1)=cx;
mu(2)=cy;
