%this function is used to make a stary image
clc;
clear;
imageSize=32;
maxGray=220;
ratio=0.95;
prepath='D:\mser\data';
fid=fopen([prepath '\realCenter.txt'],'w+');
image=zeros(imageSize,imageSize);
staryNumber=1;
for i=1:1:staryNumber

        x=6.8;
        y=7.7;
        gray=(1-rand(1)*(1-ratio))*maxGray;
        disp(gray);
        %dx=rand(1)+2;
        dx=2.3;
        dy=dx;
        pointRadius=10;
        fprintf(fid,'%d %d \t %d %d \r\n',x,y,dx,dy);
        image=makeAStart(image,gray,x,y,dx,dy,pointRadius,imageSize);
      
end
staryNumber=1
for i=1:1:staryNumber
        x=11.8;
        y=11.5;
        gray=(1-rand(1)*(1-ratio))*maxGray;
        disp(gray);
        dx=rand(1)+1.5;
        dy=dx;
        pointRadius=10;
        fprintf(fid,'%d %d \t %d %d \r\n',x,y,dx,dy);
        image=makeAStart(image,gray,x,y,dx,dy,pointRadius,imageSize);
end


staryNumber=1
for i=1:1:staryNumber

        x=20.1;
        y=24.3;
        gray=(1-rand(1)*(1-ratio))*maxGray;
        disp(gray);
        d=rand(1)+1;
        pointRadius=5;
        %用windows打开的txt。无法解析\n。只能搞成\r\n才是会
        fprintf(fid,'%d %d \t %d %d \r\n',x,y,d,d);
        image=makeAStart(image,gray,x,y,d,d,pointRadius,imageSize);
      
end
fclose(fid);

%噪点
image(26,6)=90;
% image(26,5)=250;
% image(25,5)=250;
% image(25,6)=250;

%background
image=addnoise(image,20,5);

write_infile([prepath '\star_noise_doubleStar.txt'],uint16(image));

imshow(uint8(image));
image=uint8(image);
imwrite(image,[prepath '\star_noise_doubleStar.bmp']);