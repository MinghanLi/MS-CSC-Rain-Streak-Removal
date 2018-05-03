clc; clear; close all
currentFolder = pwd; addpath(genpath(currentFolder))

load park
[h,w,~,n] = size(input);
for i=1:n
    inputY(:,:,:,i) = rgb2ycbcr(input(:,:,:,i));
end
D = reshape(inputY(:,:,1,:),[h,w,n]); D=double(D)/255;
%% initial mask param
param.lambda = 1;                                                           % lambda controls the thresholding, smaller lambda means more mask
param.weight = 1; 
param.rho =0.6; par.r=1;
par.f_size =[13 9 5];     
par.MaxIter = 1; par.method =1;                                             % 1 denotes WLRA; 2 denotes Efficient MC
par.b = 1e-2*zeros(size(par.f_size)); 
par.Mask = 1;                                                               % 0 denotes no moving objects  
par.Flam = 10; 

tic;[B, Rain, F, RainS, Filters, Mask] = CSCderain(D, param, par); toc;

DeRain = ~Mask.*B + Mask.*F;
inputY(:,:,1,:) = uint8(DeRain*255);
for i=1:n
    OutDeRain(:,:,:,i) = ycbcr2rgb(inputY(:,:,:,i));
end
I = [input(:,:,:,1:n),OutDeRain]; implay(I)
R = [Rain,RainS(:,:,:,1),RainS(:,:,:,2),RainS(:,:,:,3)]*30; implay(R)
figure;imshow(ShowFilters(Filters, 255, 0)./255);
