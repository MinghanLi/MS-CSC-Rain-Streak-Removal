function PSNR = mypsnr(distImg, origImg)
        origImg = double(origImg);
        distImg = double(distImg);

        SE= (origImg - distImg).^2;
        MSE = mean(SE(:));

        if(MSE > 0)
            PSNR = 20*log10(max(origImg(:)))-10*log10(MSE);
        else
            PSNR = 99;
        end 
end
      