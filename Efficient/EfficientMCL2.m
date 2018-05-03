%%%%% Weighted L2 norm factorization %%%%%
function [OutU,OutV,t] = EfficientMCL2(Matrix, W, InU,InV, MaxIt, tol)
%%%%% Matrix - data matrix for factorization (d*n) %%%%%
%%%%% W - weight matrix (d*n) %%%%%
%%%%% InU,InV - initialization of U (d*k dimensional) and V (n*k dimensional)     %%%%%
%%%%% MaxIt - Maximal iteration number      %%%%%
%%%%% tol - Tolerance for RobustL1 algorithm       %%%%%
 
[d n] = size(Matrix);
[k] = size(InU,2);
OutU = InU;
OutV = InV;
t = sum(sum((W.*Matrix).^2));
 
for  i = 1:MaxIt
     ind = randperm(k);
%    ind = 1:k;
    for j = ind
        TX = Matrix - OutU*OutV' + OutU(:,j)*OutV(:,j)';
        u = InU(:,j);
        [Tempv,indv] = optimMCL2(TX,W,u);    
        OutV(indv,j)=Tempv(indv);
        [Tempu,indu]= optimMCL2(TX',W',OutV(:,j));    
        OutU(indu,j) =Tempu(indu);
    end
    t = [t sum(sum(((W.*(Matrix-OutU*OutV')).^2)))];
    sum(sum(((W.*(Matrix-OutU*OutV')).^2)))/sum(sum(((W.*(Matrix)).^2)));
    
    if norm(InU - OutU) < tol
        break;
    else
        InU = OutU;
    end
end
 
Nu = sqrt(sum(OutU.^2))';
Nv = sqrt(sum(OutV.^2))';
No = diag(Nu.*Nv);
OutU = OutU*diag(1./Nu)*sqrt(No);
OutV = OutV*diag(1./Nv)*sqrt(No);

