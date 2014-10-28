%this function is used to make a stary image
%prepath='D:\mser\data';
% x=6.5;
%y=5.5
%    x=10.8;
       % y=10.5;
function image=smallstary_two(prepath,name,x1,y1,x2,y2,noiseMean,noiseVariance)
imageSize=16;
maxGray=240;
ratio=0.95;

fid=fopen([prepath '\realCenter.txt'],'w+');
image=zeros(imageSize,imageSize);
% staryNumber=1
% for i=1:1:staryNumber
% 
%         x=int32(rand(1)*imageSize);
%         y=int32(rand(1)*imageSize);
%         gray=(1-rand(1)*(1-ratio))*255;
%         d=rand(1);
%         pointRadius=5;
%       image=makeAStart(image,gray,x,y,d,pointRadius,imageSize);
% end
staryNumber=1;
for i=1:1:staryNumber

       
        x=x1;
       
        y=y1;
        gray=(1-0.8*(1-ratio))*maxGray;
        %gray=(1-rand(1)*(1-ratio))*maxGray;
        disp(gray);
        dx=1.5;
        dy=dx;
        pointRadius=10;
        %用windows打开的txt。无法解析\n。只能搞成\r\n才是会
         fprintf(fid,'%d %d \t %d %d \r\n',x,y,dx,dy);
        image=makeAStart(image,gray,x,y,dx,dy,pointRadius,imageSize);
      
end
staryNumber=1;
for i=1:1:staryNumber
        x=x2;
        y=y2;
        gray=(1-0.8*(1-ratio))*maxGray;
           %gray=(1-rand(1)*(1-ratio))*maxGray;
        disp(gray);
        dx=2;
        dy=dx;
        pointRadius=10;
         fprintf(fid,'%d %d \t %d %d \r\n',x,y,dx,dy);
        image=makeAStart(image,gray,x,y,dx,dy,pointRadius,imageSize);
end
fclose(fid);
image=addnoise(image,noiseMean,noiseVariance);
write_infile([prepath '\' name '.txt'],uint16(image));
imshow(uint8(image));
imwrite(uint8(image),[prepath '\two\' int2str(noiseMean) name '.bmp']);