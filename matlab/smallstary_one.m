%this function is used to make a stary image
% %prepath='D:\mser\data';
%  x=6.8;
% y=7.7;
function image=smallstary_one(prepath,name,x1,y1,NoiseMean,NoiseVariance)
imageSize=16;
maxGray=222;
ratio=0.95;

fid=fopen([prepath '\realCenter.txt'],'w+');
image=zeros(imageSize,imageSize);
staryNumber=1;
for i=1:1:staryNumber
       x=x1;
       y=y1;
        gray=(1-0.8*(1-ratio))*maxGray;
%         disp(gray);
        %dx=rand(1)+2;
        dx=2.3;
        dy=dx;
        pointRadius=10;
%         fprintf(fid,'%d %d \t %d %d \r\n',x,y,dx,dy);
        fclose(fid);
        image=makeAStart(image,gray,x,y,dx,dy,pointRadius,imageSize);
end

%background
image=addnoise(image,NoiseMean,NoiseVariance);
write_infile([prepath '\' name '.txt'],uint16(image));
imshow(uint8(image));

imwrite(uint8(image),[prepath '\' 'one' '\' int2str(NoiseMean) name '.bmp']);