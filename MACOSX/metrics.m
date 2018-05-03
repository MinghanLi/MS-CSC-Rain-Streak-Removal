%%range is [0 255]
function metric = metrics(groundtruth,derain)
metric = zeros(1,5);
n = size(derain,4);
rgb=size(derain,3);
% PSNR
for i = 1:rgb
    for j =1:n
        psnr(i,j) = mypsnr(derain(:,:,i,j),groundtruth(:,:,i,j));
    end
end
metric(1) = mean(psnr(:));

% VIF
for j = 1:rgb
    for i =1:n
        vif(j,i) = vifvec(double(groundtruth(:,:,j,i)),double(derain(:,:,j,i)));
    end
end
metric(2) = mean(vif(:));

% SSIM
% for j=1:rgb
%     for i = 1:n
%         ssim1(j,i) = myssim(derain(:,:,j,i),groundtruth(:,:,j,i));
%     end
% end
% metric(3) = mean(ssim1(:));

for i=1:n
    ssim2(i) = ssim(derain(:,:,:,i),groundtruth(:,:,:,i));
end
metric(5) = mean(ssim2);

% FSIM
for j=1:rgb
    for i = 1:n
        [~,fsim(i)] = FeatureSIM(groundtruth(:,:,:,i),derain(:,:,:,i));
    end
end
metric(3) = mean(fsim(:));

% UQI
for j=1:rgb
    for i = 1:n
        uqi(j,i) = img_qi(groundtruth(:,:,j,i),derain(:,:,j,i));
    end
end
metric(4) = mean(uqi(:));
