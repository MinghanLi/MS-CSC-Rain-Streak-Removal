function []=togetGif(D, SaveRoad, Correction, speed)
%”√¿¥ª≠≥ˆGIFÕº
sizeD  = size(D);
dimD   = length(sizeD);
nframe = sizeD(end);
figure
if nargin < 3
    Correction =0;
end
if nargin < 4
    speed      =0.05;
end
% if Correction
%     minD = min(D(:));
%     maxD = max(D(:));
%     D = (D-minD)/(maxD-minD);
% end

if dimD==3
    
    for i = 1:nframe
        
        imshow(abs(D(:,:,i)))  ;
        
        frame=getframe;
        im=frame2im(frame);
        [I,map]=rgb2ind(im,256);
        pause(1/24);
        if 1==i
            imwrite(I,map,[SaveRoad '.gif'],'gif','Loopcount',1,'DelayTime',speed);
        else
            imwrite(I,map,[SaveRoad '.gif'],'gif','WriteMode','append','DelayTime',speed);
        end
        
    end
    
else
    for i = 1:nframe
        
        imshow(abs(D(:,:,:,i)))  ;       
        frame=getframe;
        im=frame2im(frame);
        [I,map]=rgb2ind(im,256);
        pause(1/24);
        if 1==i
            imwrite(I,map,[SaveRoad '.gif'],'gif','Loopcount',inf,'DelayTime',speed);
        else
            imwrite(I,map,[SaveRoad '.gif'],'gif','WriteMode','append','DelayTime',speed);
        end
        
    end
    
end

end