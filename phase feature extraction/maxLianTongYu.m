%function [img]=maxLianTongYu(I) 求图像中最大的连通域
%输入：I    输入图像
%输出：img  仅包含最大连通域的图像 
function [CM]=maxLianTongYu(I)
if length(size(I))>2
    I=rgb2gray(I);
end
if ~islogical(I)
    imBw=im2bw(I);                        %转换为二值化图像
else
    imBw=I;
end
imBw=im2bw(I);                        %转换为二值化图像
imLabel=bwlabel(imBw);                %对各连通域进行标记
stats=regionprops(imLabel,'Area');    %求各连通域的大小
area=cat(1,stats.Area);
index=find(area==max(area));        %求最大连通域的索引
CM=ismember(imLabel,index);          %获取最大连通域图像