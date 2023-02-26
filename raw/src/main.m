clear
close all
clc
addpath('../utils');

filename = '0009.jpg';

I = imread(strcat('..\images\degraded\' ,filename));
I = im2gray(I);
I = im2double(I);
I = norm01(I);

% 清晰图 计算PSNR等
I_clear = imread(strcat('..\images\clear\',filename));
I_clear = im2gray(I_clear);
I_clear = im2double(I_clear);
I_clear = norm01(I_clear);

% 
% try
%    gcp('nocreate').Connected;
% catch exception
%    disp("The first run will launch the Matlab Parallel Pool, may takes few time.");
% end


tic;
I_fix = IR_correction(I);
runtime = toc;

deg_psnr = psnr(norm01(I), norm01(I_clear));
deg_ssim = ssim(norm01(I), norm01(I_clear));

res_psnr = psnr(norm01(I_fix), norm01(I_clear));
res_ssim = ssim(norm01(I_fix), norm01(I_clear));

%% plot result
disp(strcat("PSNR_degraded: ", string(deg_psnr)));
disp(strcat("SSIM_degraded: ", string(deg_ssim)));
disp("---------");
disp(strcat("PSNR_correction: ", string(res_psnr)));
disp(strcat("SSIM_correction: ", string(res_ssim)));
disp("---------");
subplot(2,2,1)
imshow(I);
title("Degraded image");
subplot(2,2,2)
imshow(I_fix)
title("Correction image");
subplot(2,2,3)
imshow(I_clear)
title("Clear image");