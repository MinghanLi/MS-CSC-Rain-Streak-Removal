function [ Mat ] = Normalize( Mat )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
dim = size(Mat,1);
den = sqrt(sum(Mat.^2))+eps;
Mat = Mat./repmat(den,dim,1);
end

