function [B, Rain,F, RainS, Filters, Mask] = CSCderain(X_Fold, param, par)  
    if (~isfield(param,'tol'))
        tol = 1.0e-5;
    else
        tol = param.tol;
    end

    if (~isfield(param,'lambda'))
        param.lambda = 5;
    end      

    if (~isfield(param,'sigma'))
        param.sigma = [];
    end      

    if (~isfield(param,'weight'))
        param.weight = 1;
    end
    
    % initialize data 
    X_Fold = padarray(X_Fold,[1 1 0]*(max(par.f_size)-1)/2,'symmetric','both');
    for a = 1:4
        X_Fold = edgetaper(X_Fold,fspecial('gaussian',6,1));
    end
    X = Unfold(X_Fold,size(X_Fold),3)';
    T = zeros(size(X)); PreL = zeros(size(X));
    Flam = par.Flam; 
    % rpca for initialize background B
    [B, U0,V0] = inexact_alm_rpca(X);
    B = gather(B); U0 = gather(U0); V0 = gather(V0);    
    r=par.r; rho = 1; FInd = FilterInd(par.f_size);
    
    % graph cuts initialization
    % GCO toolbox is called
    Omega = true(size(X_Fold));            % background support  
    ObjArea = sum(~Omega(:));
    minObjArea = numel(X(:,1))/1e4;        % minimum number of outliers
    Psigma = param.sigma; Plambda = param.lambda; Prho = param.rho;
    Pbeta = 0.5*(std(X(:,1)))^2;           % Start from a big value
    minbeta = 0.5*(3*std(X(:,1))/20)^2;    % lower bound: suppose SNR <= 20
    
    hMRF = GCO_Create(numel(X_Fold),2);
    GCO_SetSmoothCost( hMRF, [0 1;1 0] );  % Smooth Cost(3DTV)
    AdjMatrix = getAdj(size(X_Fold),param.weight);
    amplify = 10 * Plambda;
    GCO_SetNeighbors( hMRF, amplify * AdjMatrix );
    energy_cut = 0; energy_old = inf; converged = false;
     
    iter = 1; maxloopiter = 2;
    while ~converged && iter <= par.MaxIter   
       %% update Omega paramcters
        % estimate sigma
        disp([' iter = ' num2str(iter) '; update H ...']);
        loopiter = 1;
        while loopiter <= maxloopiter
        E = X-B;
        if par.Mask == 0
            Omega = ones(size(X));
        else
        if isempty(Psigma)
            sigma_old = Psigma;
            residue = sort(E(Omega(:)));
            truncate = 0.005;
            idx1 = round(truncate*length(residue))+1;
            idx2 = round((1-truncate)*length(residue));
            Psigma = std(residue(idx1:idx2));
            if abs(sigma_old-Psigma)/abs(sigma_old) < 0.01
                Psigma = sigma_old;
            end
        end
    
        % update beta
        if ObjArea < minObjArea
            Pbeta = Pbeta/2;
        else
            Pbeta =min(max([Pbeta/2,3*(Prho*Psigma)^2,minbeta]),Pbeta);
        end
        alpha = Plambda * Pbeta; 
        
        % comment these part if there is no moving object
        if Plambda > 0
            % call GCO to run graph cuts  
            GCO_SetDataCost( hMRF, int32((amplify/alpha)*[ 0.5*(E(:)).^2, ones(numel(X_Fold),1)*Pbeta ]')); 
            GCO_Expansion(hMRF);
            Omega = ( GCO_GetLabeling(hMRF) == 1 )';
            energy_cut = energy_cut + double( GCO_ComputeEnergy(hMRF) ); 
            ObjArea = sum(Omega(:)==0);
            energy_cut = (alpha/amplify) * energy_cut;
        else
            % direct hard thresholding if no smoothness
            Omega = 0.5*E.^2 < Pbeta;
            ObjArea = sum(Omega(:)==0);   
            energy_cut =  0.5*norm(X-B-E)^2+Pbeta*ObjArea;
        end
        Omega = reshape(Omega,size(X));
        end
        
        %% initial Q
        if iter==1 
           Q = Omega.*E;
        end
        
        %% Upadate U,V
        if par.method==1
            B = wlra(Omega, X-Q, B, r, 50);
        else
            [U0,V0] = EfficientMCL2(X-Q, Omega, U0,V0, 3, 1e-7);
            B = U0*V0';
        end
        difference = mean((PreL(:)-B(:)).^2);
        loopiter = loopiter+1;
        end  
        PreL = B; disp([ ' difference = ' num2str(difference)]);
        
        %% Initial Filters
        if iter == 1
            Filters = InitDictNN( reshape(Q-T, size(X_Fold)), par.f_size, FInd );     % [f_size,f_size,K]  
        end
        
 
        %% Upadate feathur map 
        disp(' Update Map ...');
        [ Map, RainS ] = SynthesisSC( reshape(Q-T,size(X_Fold)), Filters, par, FInd ); % Z[(h+f_size-1),(w+f_size-1),n,K]; RainS [h,w,n,K]  
        Rain = Unfold(sum(RainS,4),size(X_Fold),3)';
        
        
        %% Update Foreground
        disp(' Update Foreground ...');
        Temp = reshape(~Omega.*(X-Q),size(X_Fold));
        for i=1:size(X,2)
           % F_Fold(:,:,i) = FAD(Temp(:,:,i), Flam, 10, 50, [1e-4,
           % 1e-4],0);  % Flam = 1e-4
            F_Fold(:,:,i) = TVL1denoise(Temp(:,:,i), Flam, 50);
        end
        F = reshape(F_Fold,size(X));
        
        %% update Q
        TempQ = Rain+T;
        QL = (rho*TempQ+2*(X-B))/(2+rho);
        QF = (rho*TempQ+2*(X-F))/(2+rho);
        Q(Omega==1) = QL(Omega==1); Q(Omega==0) = QF(Omega==0);
        rho = rho*1.05;
        
        %% update T
        T = T+Rain-Q;
        
        %% update Filters
        disp(' Update filters ...');
        Filters = UpdateFilter(reshape(Q-T,size(X_Fold)), single(Map), Filters, par.f_size, FInd, 30);
        
        %% updata model param
        Error = Omega.*(X-B-Q);
        par.b = mean(abs(reshape(Map,numel(X),length(par.f_size))));
        energy = energy_cut + sum(Error(:).^2);
        if ObjArea > minObjArea && abs(energy_old-energy/energy) < tol; break; end
        energy_old = energy;
        if (difference<1e-30); iter = par.MaxIter+1; end
        iter = iter+1;
    end
    maxf_size =max( par.f_size);
    B = EdgaCut(reshape(B,size(X_Fold)),maxf_size);
    Rain = EdgaCut(reshape(Rain,size(X_Fold)),maxf_size);
    F = EdgaCut(reshape(F,size(X_Fold)),maxf_size);
    Map = EdgaCut(Map,maxf_size);                           % [h,w,rgb,n,k]
    Mask = EdgaCut(reshape(~Omega,size(X_Fold)),maxf_size); 
    RainS = EdgaCut(RainS,maxf_size);
    Q = EdgaCut(reshape(Q,size(X_Fold)),maxf_size);
end

%% function to get the adjacent matirx of the graph
function W = getAdj(sizeData,weight)
numSites = prod(sizeData);
id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 1+1:numSites+1,...
1+sizeData(1):numSites+sizeData(1),...
1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
%value = ones(1,3*numSites);
value = [weight*ones(1,2*numSites),1*ones(1,numSites)];
W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);
end

function Ind = FilterInd(f_size)
K = length(f_size); maxf_size = max(f_size);
ind = (maxf_size-f_size)/2;
IndMax = zeros(maxf_size,maxf_size,K);
for k=1:K
    IndMax(1+ind(k):end-ind(k),1+ind(k):end-ind(k),k) = 1;
end
Ind = IndMax==1; 
end

function  Filters  = InitDictNN( R, f_size, FInd )
% pca for initialize D
maxf_size = max(f_size); K = length(f_size);
Patches = Video2Patch( R, maxf_size );
Patches = Patches- repmat(mean(Patches,2),1,size(Patches,2)); 
temp = Patches*Patches';
[U, ~, ~] = svd(temp);
D = max(U(:,2:K+1),0); 
Temp = zeros(size(D)); Temp(FInd) = D(FInd);
Filters = reshape(Normalize(Temp),[maxf_size,maxf_size,K]);
end

function Filters = UpdateFilter( X, Z, oriF, f_size, FInd, MaxIter )
% UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
[maxf_size,~,K] = size(oriF);
X = EdgaCut(X,maxf_size);      
f_vec = oriF(FInd)';
Zpat = Video2Patch(Z, f_size);
X_vec = X(:)';

iter=1; step = 1/100;
while(iter < MaxIter)
    Grad  = Zpat*(f_vec*Zpat - X_vec)';    
    f_vec = max(f_vec-step*Grad',1e-5);
    cowInd = 0;
    for k=1:K
        f_vec(cowInd+1:cowInd+f_size(k)^2) = Normalize(f_vec(cowInd+1:cowInd+f_size(k)^2));
        cowInd = cowInd+f_size(k)^2;
    end
    iter = iter+1;
    if  norm(Grad)<1e-5
        break;
    end
end
Filters = zeros(size(oriF));
Filters(FInd) = f_vec;
end

function [ M, Rain ] = SynthesisSC( X, Filters, par, FInd )
% UNTITLED2 Summary of this function goes here
% Detailed explanation goes here
[h,w,n] = size(X); mu=par.b; f_size=par.f_size;
for k=1:size(Filters,3)
    Fk = Filters(:,:,k);Temp = reshape(Fk(FInd(:,:,k)),[f_size(k),f_size(k)]);
    FFilters(:,:,:,k) = psf2otf(rot90(Temp,2),[h,w,n]);                     % [h+2*CutEdge,w+2*CutEdge,K]
end
[M, Rain]= CSC_ADMM_GPU( fft2(X), single(FFilters), mu, 200);
end


