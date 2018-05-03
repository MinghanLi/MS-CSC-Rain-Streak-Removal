%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code implements the algorithm proposed in                        %% 
%% A Highly Parallel Algorithm for Isotropic Total Variation Models by   %%
%% Jie Wang, Qingyang Li, Sen Yang, Wei Fan, Peter Wonka and Jieping Ye. %%
%% The algorithm aims to solve the isotropic ROF model:                  %%
%% min_X 0.5|X - Y|_F^2 + lambda|X|_TV,                                  %%
%% where the TV norm is the summation of the l2 norm of the gradient at  %%
%% point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%% setup the model parameters
lam = 0.35; % the regularization parameter

%% setup the parameters of FAD
MAX_ITER = 100; % the maximum iteration number

TOL_SETTING = [1e-4,1e-4]; % tolerance for the stopping condition

DISPLAY_SETTING = 1; % display the errors in each iteration; otherwise, set it to 0

gamma = 10; % the dual accent step length

%% generate the synthetic image and add noise
rn = 0.2;   % std of the Gaussian noise
%I = im2double(imread('Brain.png')); % if the image is gray scale
I = im2double(rgb2gray(imread('Lenna.png'))); % if the image has RGB colors.
Im = I + rn*randn(size(I)); % the noised image

%% start 

t = [];

fprintf('FAD is running...\n');
tstart = tic;
[X, fun, iterNum] = FAD(Im,lam,gamma,MAX_ITER,TOL_SETTING,DISPLAY_SETTING);
tloc=toc(tstart);
t=[t tloc];

fprintf('\npFAD is running...\n');
tstart = tic;
[X1, fun1, iterNum1] = pFAD_OMP(Im,lam,gamma,MAX_ITER,TOL_SETTING,DISPLAY_SETTING);
tloc=toc(tstart);
t = [t tloc];


%% plot
figure(1);
subplot(2,2,1); imshow(I); title('Original');
subplot(2,2,2); imshow(Im); title('noisy');
subplot(2,2,3); imshow(X); title('denoised by FAD');
subplot(2,2,4); imshow(X1); title('denoised by pFAD\_OMP');






