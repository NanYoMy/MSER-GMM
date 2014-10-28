%用于显示高度图
prepath='D:\mser\data';
center=load([prepath '\star.txt']);
[h,w]=size(center);
surf([1:h],[1:w],center);
colormap([1 1 0;1 1 1]);
grid off;
axis off;