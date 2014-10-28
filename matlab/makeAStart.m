function image=makeAStart(img,gray,x,y,dx,dy,pointRadius,imageSize)


% imageSize=1024;
% pointRadius=3;
% gray=255;
% x=50;
% y=50;
% d=1;
% img=zeros(imageSize,imageSize);

%几何中心，为质心。

x=double(x);
y=double(y);
dx=double(dx);
dy=double(dy);
gray=double(gray);
pointRadius=double(pointRadius);
startX=int32(x-pointRadius)-1;
startY=int32(y-pointRadius)-1;
endX=int32(x+pointRadius)+1;
endY=int32(y+pointRadius)+1;
%testPoint=zeros(pointRadius*2+1,pointRadius*2+1);
% testPoint=zeros(imageSize,imageSize);
%search the point
for i=startX:1:endX
    for j=startY:1:endY
        %check the bounder
        if i<1|i>imageSize|j<1|j>imageSize
            continue;
        end
        
        %the stat must be a circle
        distance=(i-x-0.5)*(i-x-0.5)+(j-y-0.5)*(j-y-0.5);
        if distance>(pointRadius*pointRadius)
            continue;
        end
        
%         if (i==9 &j==10)|(i==10 & j==10)|(i==10 & j==9)
%             display(i);
%         end
        
        %guassion
  
        a=-double((double(i)-x-0.5)*(double(i)-x-0.5)+(double(j)-y-0.5)*(double(j)-y-0.5));
        temp=a/2/dx/dy;
        temp=double(temp);
        
        img(i,j)=((exp(temp)/2/pi/dx/dy)*gray*2*pi*dx*dy+img(i,j));
%         disp(img(i,j));
%         img(i,j)=(exp(temp)/2/pi/d/d)*gray*2*pi*d*d+(rand(1)*2-1)*0.05*255;
%         testPoint(i-x+pointRadius+1,j-y+1+pointRadius)=(exp(temp)/2/pi/d/d)*gray*2*pi*d*d;
    end
end
image=img;
% imShow(img);
%X=1:1:pointRadius*2+1;
%Y=1:1:pointRadius*2+1;
%imshow(uint8(testPoint))
%surf(X,Y,testPoint)