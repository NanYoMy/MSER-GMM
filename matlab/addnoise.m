function noiseimg = addnoise(img, mean, sigma)
noiseimg = img;
[h w] = size(noiseimg);

for i = 1:h
    for j = 1:w
        n = mean + sigma.*randn;
        noiseimg(i,j) = noiseimg(i,j) + n;
    end
end

return