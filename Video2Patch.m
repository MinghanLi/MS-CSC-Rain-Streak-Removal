 function  Patches  =  Video2Patch( D,f_size )
maxf_size = max(f_size); r=0;
TotalPatNum = (size(D,1)-maxf_size+1)*(size(D,2)-maxf_size+1)*size(D,3);                  %Total Patch Number in the image
Patches   =   zeros(sum(f_size.^2), TotalPatNum,'single');                      %Current Patches 
%D = gpuArray(D); 
ind = (maxf_size-f_size)/2;
for k = 1:length(f_size)
    PatSize = f_size(k);
    for j = 1:PatSize
        for i  = 1:PatSize
            r=r+1;
            temp =  D(i+ind(k):end-ind(k)-PatSize+i,j+ind(k):end-ind(k)-PatSize+j,:,k);  % 按列展开[PatSize*PatSize, TotalPatNum]     
            Patches(r,:)=  temp(:)';
        end
    end 
end
%Patches = gather(Patches);
 