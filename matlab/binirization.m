%??????

%prepath
function biImage= binirization(prepath,name,image)
maxT=max(max(image));
minT=min(min(image));
[h,w]=size(image);
%ÅÅÐò·Ö¸î
threadCaculateArray=reshape(image,1,h*w);
threadCaculateArray=sort(threadCaculateArray);
threadCaculateArray=threadCaculateArray((h*w*0.2):(h*w*0.7));
mu=mean(threadCaculateArray);
sigma=std(threadCaculateArray);
t=mu+3*sigma;

tmp=image;

for i=1:size(tmp,1)
    for j=1:size(tmp,2)

        if tmp(i,j)>t
            tmp(i,j)=1;
        else
            tmp(i,j)=0;
        end
    end
end
biImage=tmp;
imshow(tmp);

write_infile([prepath '\' name  '.txt'],uint16(biImage));
imwrite(biImage,[prepath '\' name '.bmp']);



