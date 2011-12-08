function [k,I] = coarse_kernel_est(Ish_x,Ish_y,Im_x,Im_y,Im,ksize,lambda,ganma)
%% This function implement coarse-to-fine kernel estimation
% E(k)=||diff_Ish(*)k-diff_Im||^2+ganma*||k||^2
% E(I)=||I(*)k-Im||^2+lambda*||diff_I-diff_Ish||^2
% (*):represent the convolution
% diff_*:represent the first-order derivative

[xim,yim]=size(Im);
F1 = [1,-1];
F2 = [1;-1];

FFtIsh_x = fft2(Ish_x);
FFtIsh_y = fft2(Ish_y);
FFtIm_x = fft2(Im_x);
FFtIm_y = fft2(Im_y);
FFtF1 = psf2otf(F1,[xim,yim]);
FFtF2 = psf2otf(F2,[xim,yim]);
FFtIm = fft2(Im);

Nomin = conj(FFtIsh_x).*FFtIm_x+conj(FFtIsh_y).*FFtIm_y;
Denom = conj(FFtIsh_x).*FFtIsh_x+conj(FFtIsh_y).*FFtIsh_y+ganma;
FFtk = (Nomin./Denom);
k = otf2psf(FFtk,[ksize,ksize]);
k(find(k <= 0)) = 0;
ksum = sum(k(:));
k = k./ksum;

FFtk = psf2otf(k,[xim,yim]);
Nomin = conj(FFtk).*FFtIm+lambda*(conj(FFtF1).*FFtIsh_x+conj(FFtF2).*FFtIsh_y);
Denom = conj(FFtk).*FFtk+lambda*(conj(FFtF1).*FFtF1+conj(FFtF2).*FFtF2);
I = ifft2(Nomin./Denom);