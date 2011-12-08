function yout = multi_deriv_deconv(yin, k, lambda)
%% This fucntion implement for Multi-Derivative-Deconvolution
weight=2;
weit=[weight,weight/2,weight/2,weight/4,weight/4,weight/4,weight/4];

F1=[1,-1];
F2=[1;-1];
F3=[1,-2,1];
F4=[1;-2;1];
F5=[1,-2;0,1];
F6=[1,0;-2,1];


% initialize with input or passed in initialization
yout = yin; 

% make sure k is a odd-sized
if ((mod(size(k, 1), 2) ~= 1) | (mod(size(k, 2), 2) ~= 1))
  fprintf('Error - blur kernel k must be odd-sized.\n');
  return;
end;
ks = floor((size(k, 1)-1)/2);

% compute constant quantities
% see Eqn. (3) of paper
[Nomin1, Denom1, Denom2] = computeDenominator(yin, k,weit,F1,F2,F3,F4,F5,F6);

gamma =  lambda;
Denom = Denom1 + gamma*Denom2;
Fyout = Nomin1./Denom; 
yout = real(ifft2(Fyout));
      
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nomin, Denom, Denom2] = computeDenominator(y, k,weit,F1,F2,F3,F4,F5,F6)
%
%
% Inputs: 
%  y: blurry and noisy input
%  k: convolution kernel  
% 
% Outputs:
%      Nomin1  -- F(K)'*F(y)
%      Denom1  -- |F(K)|.^2
%      Denom2  -- |F(D^1)|.^2 + |F(D^2)|.^2
%

sizey = size(y);
otfk  = psf2otf(k, sizey); 

FF1=psf2otf(F1,sizey);
FF2=psf2otf(F2,sizey);
FF3=psf2otf(F3,sizey);
FF4=psf2otf(F4,sizey);
FF5=psf2otf(F5,sizey);
FF6=psf2otf(F6,sizey);
Y=fft2(y);
% if higher-order filters are used, they must be added here too
Nomin = weit(1)*conj(otfk).*Y;
Nomin = Nomin+weit(2)*conj(FF1).*conj(otfk).*FF1.*Y; 
Nomin = Nomin+weit(3)*conj(FF2).*conj(otfk).*FF2.*Y; 
Nomin = Nomin+weit(4)*conj(FF3).*conj(otfk).*FF3.*Y; 
Nomin = Nomin+weit(5)*conj(FF4).*conj(otfk).*FF4.*Y; 
Nomin = Nomin+weit(6)*conj(FF5).*conj(otfk).*FF5.*Y; 
Nomin = Nomin+weit(7)*conj(FF6).*conj(otfk).*FF6.*Y; 

Denom = weit(1)*conj(otfk).*otfk;
Denom = Denom+weit(2)*conj(FF1).*conj(otfk).*otfk.*FF1; 
Denom = Denom+weit(3)*conj(FF2).*conj(otfk).*otfk.*FF2; 
Denom = Denom+weit(4)*conj(FF3).*conj(otfk).*otfk.*FF3; 
Denom = Denom+weit(5)*conj(FF4).*conj(otfk).*otfk.*FF4; 
Denom = Denom+weit(6)*conj(FF5).*conj(otfk).*otfk.*FF5; 
Denom = Denom+weit(7)*conj(FF6).*conj(otfk).*otfk.*FF6; 

Denom2 = conj(FF1).*FF1 + conj(FF2).*FF2;
