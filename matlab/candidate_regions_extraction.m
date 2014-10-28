%candadidate region extractin

%
%prepath='D:\mser\data';
function candidate=candidate_regions_extraction(prepath,name,originalImage,biImage)

image=originalImage;
imagePattern=biImage;
[h,w]=size(image);
candidate=zeros(h,w);
for i=1:h
    for j=1:w
    
        if imagePattern(i,j)==1
            candidate(i,j)=image(i,j);
        else
            candidate(i,j)=0;
        end
        
    end
end
write_infile([prepath '\' name '.txt'],uint16(candidate));

imshow(uint8(candidate));

imwrite(uint8(candidate),[prepath '\' name '.bmp']);