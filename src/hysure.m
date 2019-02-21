function [SRI_hat,info] = hysure(HSI,MSI, opt)
%UNTITLED23 Summary of this function goes here
%   Detailed explanation goes here

% hsize_h = 10; hsize_w = 10;
% shift = 1;
hsize_h = 10; hsize_w = 10;
shift = 3; blur_center = 4; lambda_R = 1e1; lambda_B = 1e1;
basis_type = 'VCA'; lambda_phi = 5e-4; lambda_m = 1e0;

d1 = ceil(size(MSI,1)/size(HSI,1));

if size(MSI,3)==4
    Km = 4;
    spectral_distrib = linspace(430,860,size(HSI,3));
    bands = [430 545; 466 620; 590 710; 715 918];
    Pm = zeros(Km,size(HSI,3)); intersection = {};
    for k=1:Km
        ind = find(spectral_distrib >= bands(k,1) & spectral_distrib <= bands(k,2));
        intersection{k} = ind;
    end
    contiguous = intersection;
elseif size(MSI,3)==6
    Km = 6;
    spectral_distrib = linspace(400,2500,size(HSI,3)); 
    bands = [450 520; 520 600; 630 690; 760 900; 1550 1770; 2080 2350];
    Pm = zeros(Km,size(HSI,3)); intersection = {};
    for k=1:Km
        ind = find(spectral_distrib >= bands(k,1) & spectral_distrib <= bands(k,2));
        intersection{k} = ind;
    end
    contiguous = intersection;
elseif size(MSI,3)==28
    d = size(HSI,3)/size(MSI,3);
    for k=1:28
        intersection{k} = (k-1)*d+1:k*d;
    end
    contiguous = intersection;
elseif size(MSI,3)==1
    intersection{1} = 1:size(HSI,3);
    contiguous = intersection;
end

%[V, R_est, B_est] = sen_resp_est(HSI, MSI, d1, intersection, contiguous, opt.p, lambda_R, lambda_B, hsize_h, hsize_w, shift, blur_center);
%tic; SRI_hat = data_fusion(HSI, MSI, d1, R_est, B_est, opt.p, basis_type, lambda_phi, lambda_m); t = toc;

Pm = zeros(1,size(HSI,3));
spec_range = linspace(400,2500,size(HSI,3));
Pm(1,:) = 1/length(spec_range);
R_est = Pm;
B_est = zeros(size(MSI,1),size(MSI,2));
h1 = gauss_kernel(9, 1);
B_est(1:9,1:9) = h1' * h1;

%B_est = circshift(circshift(B_est,-4,2),-4,1);

tic; SRI_hat = data_fusion(HSI, MSI, d1, R_est, B_est, opt.p, basis_type, lambda_phi, lambda_m); t = toc;
%Zhat = im2mat(SRI_hat); Zhat_denoised = (V*V')*Zhat; SRI_hat = mat2im(Zhat_denoised, size(MSI,1));

Zhat = im2mat(SRI_hat); Zhat_denoised = Zhat; SRI_hat = mat2im(Zhat_denoised, size(MSI,1));
info.time = t;

end

