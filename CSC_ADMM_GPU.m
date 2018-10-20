function [M, Rain] = CSC_ADMM_GPU( X, F, lambda, MaxIter)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
max_rho = 1e5; rho=1; tol = 1e-10;
[h,w,n] = size(X); K=size(F,3);
iter = 1; Cond = 1;

lambda=gpuArray(repmat(reshape(lambda,[1,1,1,K]),[h,w,n,1]));
F = gpuArray(single(F));
for k=1:K
    Filters(:,:,:,k) = psf2otf(rot90(F(:,:,k),2),[h,w,n]);
end
FX = gpuArray(fft2(single(X)));
MU =  gpuArray.zeros(size(Filters),'single');
Y =  gpuArray.zeros(size(Filters),'single');

C_Filters = conj(Filters);            % [h,w,n,K]
FTX = C_Filters.*repmat(FX,[1,1,1,K]);  % [h,w,n,K]
FTF  = sum(C_Filters.*Filters,4);     % [h,w,n]

while(iter < MaxIter && Cond)
    FR = FTX+rho*fft2(Y-MU);            % [h,w,n,K]
    FM = (FR-repmat(sum(Filters.*FR,4)./(rho+FTF), [1,1,1,K]) .*C_Filters)./rho;  % [h,w,n,K]    
    M = real(ifft2(FM));   % [h,w,n,K]
    Y = max(M+MU-lambda/rho,0);  % [h,w,n,,K]
    MU = MU + M - Y;
    if(rho < max_rho)
        rho = rho*1.01;
    end
    if(mod(iter,10)==0)
        ConvergenceError = mean((M(:)-Y(:)).^2);
        Cond = (ConvergenceError>tol);
    end
    iter = iter+1;
end
M = double(gather(Y));
Rain = double(gather(real(ifft2(FM.*Filters))));    


