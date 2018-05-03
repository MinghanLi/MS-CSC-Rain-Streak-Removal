function newim = TVL1denoise(im, lambda, niter)
%figure; imshow(outim, []);

% TVL1denoise - TV-L1 image denoising with the primal-dual algorithm.
%
% Function denoises an image using the TV-L1 model optimized with a primal-dual
% algorithm; this works best for impulse noise. For more details, see
%   A. Mordvintsev: ROF and TV-L1 denoising with Primal-Dual algorithm,
%   http://znah.net/rof-and-tv-l1-denoising-with-primal-dual-algorithm.html
% and
%   Chambolle et al. An introduction to Total Variation for Image Analysis, 2009. <hal-00437581>
%   https://hal.archives-ouvertes.fr/hal-00437581/document
%
% Usage:
%    newim = TVL1denoise(im, lambda, niter)
%
% Arguments:
%    im          -  Image to be processed
%    lambda      -  Regularization parameter controlling the amount
%                   of denoising; smaller values imply more aggressive
%                   denoising which tends to produce more smoothed results.
%    niter       -  Number of iterations. If omitted or empty defaults
%                   to 100
%
% Returns:
%    newim       -  Denoised image with pixel values in [0 1].

% Copyright (c) 2016 Manolis Lourakis (lourakis **at** ics forth gr)
% Institute of Computer Science,
% Foundation for Research & Technology - Hellas
% Heraklion, Crete, Greece.
% http://www.ics.forth.gr/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

  if(nargin<3) || isempty(niter)
    niter=100;
  end

  L2=8.0;
  tau=0.02;
  sigma=1.0/(L2*tau);
  theta=1.0;
  lt=lambda*tau;

  [height width]=size(im); 
  unew=zeros(height, width);
  p=zeros(height, width, 2);
  d=zeros(height, width);
  ux=zeros(height, width);
  uy=zeros(height, width);

  mx=max(im(:));
  if(mx>1.0)
    nim=double(im)/double(mx); % normalize
  else
    nim=double(im); % leave intact
  end
  u=nim;

  %[p(:, :, 1), p(:, :, 2)]=imgradientxy(u, 'IntermediateDifference');
  p(:, :, 1)=u(:, [2:width, width]) - u;
  p(:, :, 2)=u([2:height, height], :) - u;

  for k=1:niter
    % projection
    % compute gradient in ux, uy
    %[ux, uy]=imgradientxy(u, 'IntermediateDifference');
    ux=u(:, [2:width, width]) - u;
    uy=u([2:height, height], :) - u;
    p=p + sigma*cat(3, ux, uy);
    % project
    normep=max(1, sqrt(p(:, :, 1).^2 + p(:, :, 2).^2)); 
    p(:, :, 1)=p(:, :, 1)./normep;
    p(:, :, 2)=p(:, :, 2)./normep;

    % shrinkage
    % compute divergence in div
    div=[p((1:height-1), :, 2); zeros(1, width)] - [zeros(1, width); p((1:height-1), :, 2)];
    div=[p(:, (1:width-1), 1)  zeros(height, 1)] - [zeros(height, 1)  p(:, (1:width-1), 1)] + div;

    %% TV-L2 model
    %unew=(u + tau*div + lt*nim)/(1+tau);

    % TV-L1 model
    v=u + tau*div;
    unew=(v-lt).*(v-nim>lt) + (v+lt).*(v-nim<-lt) + nim.*(abs(v-nim)<=lt);
    %if(v-nim>lt); unew=v-lt; elseif(v-nim<-lt) unew=v+lt; else unew=nim; end

    % extragradient step
    u=unew + theta*(unew-u);

    %% energy being minimized
    % ux=u(:, [2:width, width]) - u;
    % uy=u([2:height, height], :) - u;
    % E=sum(sqrt(ux(:).^2 + uy(:).^2)) + lambda*sum(abs(u(:) - nim(:)));
    % fprintf('Iteration %d: energy %g\n', k, E);
  end

  newim=u;
end