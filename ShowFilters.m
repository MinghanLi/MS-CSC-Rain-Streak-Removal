function [ BigImage ] = ShowFilters( filters , scalar,flag)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
filtersize = size(filters,1);
filternum = size(filters,3);
n = size(filters,4);

BigImageWidth = (2 * filtersize + 3*5)*n;
BigImageHeight = (floor(filternum/2)+1)*(filtersize+6) + 6;

BigImage = 255*ones(BigImageHeight, BigImageWidth);
for j =1:n
    CurHeight = 10;
    CurWidth  = (j-1)*2*(filtersize+5)+5;
    for i =1: filternum
    tempfilter = reshape(filters(:,:,i,j),filtersize,filtersize);
    BigImage(CurHeight+1:CurHeight+filtersize, CurWidth+1:CurWidth+filtersize) = scalar*tempfilter+128;
    if(mod(i,2)==0)
        CurWidth = (j-1)*2*(filtersize+5)+5;
        CurHeight = CurHeight+filtersize+5;
    else
        CurWidth = CurWidth + filtersize+5;
    end
    end
end
if(flag)
imshow(BigImage/255);
end
end

