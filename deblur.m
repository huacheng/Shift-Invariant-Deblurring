%---------------------------------------------------------%
%  This program implement shift-invariant Deblurring  %
%         Creative Commons(CC) Hua Cheng, 2011             %
%---------------------------------------------------------%

close all;
clear all;

I = im2double(imread('img/blurry0.tiff'));
I = I(:,:,1);
lambda = 0.001;
ganma = 10;
tau_r = 0;
tau_s = 0;
tic

ratio = [10,7,5,3,2,1.5,1];
[~,yratio] = size(ratio);
ks = [1,2,4,7,9,12,17];
Im = imresize(I,1/ratio(1));
I_latent = Im;

ind=1:yratio;
for i=ind
    i
    %% Imresize for multi-scale computing
    
    Im = imresize(I,1/ratio(i));
    [xim,yim] = size(Im);
    [xil,yil] = size(I_latent);
    minxi = min([xim,xil]);
    minyi = min([yim,yil]);
    Im = Im(1:minxi,1:minyi);
    I_latent = I_latent(1:minxi,1:minyi);
    
    
    %% Image expanding
    ksize = 2*ks(i)+1;
    sigma = 2;
    kernel_gaus = fspecial('gaussian', ksize, sigma);
    %%edgetaper to better handle circular boundary conditions
    Im = padarray(Im, [1 1]*ks(i), 'replicate', 'both');
    I_latent = padarray(I_latent, [1 1]*ks(i), 'replicate', 'both');
    for j=1:3
        Im = edgetaper(Im, kernel_gaus);
        I_latent = edgetaper(I_latent,kernel_gaus);
    end

    % Im=imnoise(Im,'gaussian',0,0.005);
    % figure;imshow(Im)


    %% Perona and Malik nonliear Deffusion
    iter=5;
    I_sh=perona_malik(I_latent,iter);
    %% time and spatial steps
    % Oshr et.al. shock filter 
    dt=0.1; h=1;
    iter=5;
    I_sh=shock_filter(I_sh,iter,dt,h);  
    
    %% fisrt order derivative for the blur image Im and the shock filter iamge Ish
    [xim,yim]=size(Im);
    Im_x = [Im(1:xim-1,:)-Im(2:xim,:);Im(xim,:)-Im(1,:)];
    Im_y = [Im(:,1:yim-1)-Im(:,2:yim),Im(:,yim)-Im(:,1)];
    Ish_x = [I_sh(1:xim-1,:)-I_sh(2:xim,:);I_sh(xim,:)-I_sh(1,:)];
    Ish_y = [I_sh(:,1:yim-1)-I_sh(:,2:yim),I_sh(:,yim)-I_sh(:,1)];

    %% M = Heaviside(rmap-tau_r)
    [M tau_r] = M_compute(Im_x,Im_y,ks,i,tau_r);
    %% H = Heaviside(M.*||Delta_Ish||_2-tau_s)
    [H tau_s] = H_compute(Ish_x,Ish_y,M,ks,i,tau_s);
    %% \downtriangular_I^s computing
    %  \downtriangular_I^s = Delta_Ish.*H
    Is_x = Ish_x.*H;
    Is_y = Ish_y.*H;
    clear Ish_x Ish_y;
    [k,I_latent] = coarse_kernel_est(Is_x,Is_y,Im_x,Im_y,Im,ksize,lambda,ganma);
    
    if i < ind(end)
        I_latent = I_latent(ks(i)+1:xim-ks(i),ks(i)+1:yim-ks(i));
        I_latent = imresize(I_latent,ratio(i)/ratio(i+1),'bicubic');        
    end
end
toc;

% % save k
% 
% %% Iterative Surpport Detection for Kernel Refinement
% % This result has failed!!
% % I have not implement the Idea of ISD successfully!
% tic
% % load k
% k0 = k;
% figure;
% imagesc(k)
% title('K')
% colormap gray
% 
% iter = 4;
% innerIter = 1;
% beta = 1;
% epslon_s = 0;
% 
% FFtIs_x = fft2(Is_x);
% FFtIs_y = fft2(Is_y);
% FFtIm_x = fft2(Im_x);
% FFtIm_y = fft2(Im_y);
% FFtAtB = conj(FFtIs_x).*FFtIm_x;
% FFtAtB = FFtAtB+conj(FFtIs_y).*FFtIm_y;
% AtB = otf2psf(FFtAtB,[ksize,ksize]);
% FFtAtA = conj(FFtIs_x).*FFtIs_x;
% FFtAtA = FFtAtA+conj(FFtIs_y).*FFtIs_y;
% 
% for i = 1:iter
% 
%     f = zeros(ksize,ksize);
%     di = max(k0(:))/(2*ksize*i);
%     k_sort = sort(k0(:),'ascend');
%     [xk,yk] = size(k_sort);
%     for j=1:xk-1
%         if ((k_sort(j+1)-k_sort(j)) > di)
%             epslon_s = k_sort(j);
%             break;
%         end
%     end
%     for idx = 1:ksize
%         for idy = 1:ksize
%             if k0(idx,idy) <= epslon_s
%                f(idx,idy) = 1;
%             end
%         end
%     end
%     
%     k1 = k0;
%     for j=1:innerIter
%         k1_sign = sign(k1);
%         f_k_sign = f.*k1_sign;
%         FFtAtA_k = FFtAtA.*psf2otf(k1,[xim,yim]);
%         AtA_k = otf2psf(FFtAtA_k,[ksize,ksize]);
%         
%         FFtAtA_k_sign = FFtAtA.*psf2otf(k1_sign,[xim,yim]);
%         AtA_k_sign = otf2psf(FFtAtA_k_sign,[ksize,ksize]);
%         
%         E_k = AtA_k-AtB+beta*f_k_sign;
%         E_kk = AtA_k_sign;
%         
%         k1 = k1-E_k./E_kk;
%         
%         k1(find(k1 <= 0)) = 0;
%         k1_sum = sum(k1(:));
%         k1 = k1./k1_sum;
%         
% %         figure;
% %         imagesc(k1)
% %         colormap gray
% %         title('k1')
% %         colorbar
%     end
%     
%     k0 = k1;
% % 
% %     figure;
% %     imagesc(1-f)
% %     colormap gray
% %     title('f')
% %  
% 
% %     figure;
% %     imagesc(k0)
% %     colormap gray
% %     title('k0')
% %     colorbar 
% 
% end
% % 
% % toc;

tic;
lambda = 0.05;
I_latent = multi_deriv_deconv(Im, k, lambda);
toc;
% Im = Im(ks(end)+1:xim-ks(end),ks(end)+1:yim-ks(end));
I_latent = I_latent(ks(end)+1:xim-ks(end),ks(end)+1:yim-ks(end));
imwrite(I_latent,'img/blurry0d.tiff');

% figure;
% imshow(Im);
% title('B');
% 
% figure;
% imshow(I_latent)
% title('I\_latent')

% figure;
% imshow(I_sh);
% title('SH')

% figure;
% imagesc(k)
% title('K')
% colormap gray



% figure;
% imshow(I_latent0)
% title('I\_latent0')









