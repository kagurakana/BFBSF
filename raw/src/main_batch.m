clear
close all
clc
addpath('../utils');

deFloader = fullfile('../images/degraded/');
clFloader = "../images/clear/";
normlized = true;

flist = dir(fullfile(deFloader,'*.*'));
addpath('../utils');

fileNames={flist.name};
fileNames=fileNames(3:end);
i=0;
F=zeros(size(fileNames,2),1);

psnrArr_M = [];

ssimArr_M = [];

gssimArr_M = [];

timerArr_M = [];

idx = 1:length(fileNames);


for ii=1:length(fileNames)

    file_idx = idx(ii);
    p = [0,0,0,0];
    filename =strcat(deFloader, fileNames{file_idx});
    ii
    I = imread(filename);
    I= im2double(im2gray(I));
    if normlized
        I = norm01(I);
    end
    filename =strcat(clFloader, fileNames{file_idx});
    I_clear = imread(filename);
    I_clear = im2double(im2gray(I_clear));
    I_clear = norm01(I_clear);

    %     [psnrArr, ssimArr,gssimArr,timerArr] = correction_with_methods(I,I_clear,fileNames{file_idx},...
    %         {@Cao_method, @biasCorrection_bipoly,@lile_method,@ShiMethod, @cheb_method,@IR_correction},...
    %         ["cao","zheng","lile","butterworth","shi","bezier"],true);
    
    [psnrArr, ssimArr,gssimArr,timerArr] = correction_with_methods(I,I_clear,fileNames{file_idx},...
        {@IR_correction},...
        ["BFBSF"],true);
    
    psnrArr_M = [psnrArr_M;psnrArr];
    ssimArr_M = [ssimArr_M;ssimArr];
    gssimArr_M = [gssimArr_M;gssimArr];
    timerArr_M = [timerArr_M;timerArr];
end


disp('-----DONE-----');
disp('DEG PSNR:  MEAN CORR PSNR: ');
disp(sum(psnrArr_M)/length(psnrArr_M));
disp('--------------');
disp('DEG SSIM:  MEAN CORR SSIM: ');
disp(sum(ssimArr_M)/length(ssimArr_M));
disp('--------------');
disp('MEAN TIME of MSP-BFBSF: ');
disp(sum(timerArr_M)/length(psnrArr_M));


function [psnr_list, ssim_list, gssim_list, time_list] = correction_with_methods(I,I_clear,filename,methods_list,methods_name_list,is_save)
m_len = length(methods_list);
psnr_list = zeros(1,m_len +1);
ssim_list = zeros(1,m_len +1);
gssim_list = zeros(1,m_len +1);
time_list = zeros(1,m_len);
if exist('I_clear','var') && class(I_clear)~="char" && class(I_clear)~="string"
    psnr_list(1) = psnr(I,I_clear);
    ssim_list(1) = ssim(I,I_clear);
    gssim_list(1) = gssim(I,I_clear);
end

for i = 2:m_len + 1
    [~,I_out] = evalc('methods_list{i-1}(I)');
    I_out = im2gray(im2double(norm01(I_out)));
    if exist('I_clear','var')  && class(I_clear)~="char" && class(I_clear)~="string"
        psnr_list(i) = psnr(I_out,I_clear);
        ssim_list(i) = ssim(I_out,I_clear);
        gssim_list(i) = gssim(I_out,I_clear);
        time_list(i) = time_res;
    end
    if exist('is_save','var')
        save_folader = './result/';
        if ~isfolder(strcat(save_folader, methods_name_list(i - 1)))
            mkdir(strcat(save_folader, methods_name_list(i - 1)));
        end
        if filename(end-3:end) == ".jpg"
            imwrite(I_out,strcat(save_folader, methods_name_list(i - 1),'/',filename),'Quality',100)
        else
            imwrite(I_out,strcat(save_folader, methods_name_list(i - 1),'/',filename))
        end
    end
end
end