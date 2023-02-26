% Progressive BFBSF+ Algorithm
function I_out = IR_correction_complex(I,ita,alpha,beta,gamma)

addpath('../utils');

count =16;   % 迭代次数

if ~exist('ita','var')
    ita = 0.2;
end

if ~exist('alpha','var')
    alpha = 0.6;
end

if ~exist('beta','var')
    beta = 0.8;
end

if ~exist('gamma','var')
    gamma = 0.3;
end

[size_I_x,~] = size(I);

down_sample_rate = 1/2;

d = 6;
C_N_I = zeros(d);

for i = 2:d
    for j = 0:i
        C_N_I(i-1,j+1) =nchoosek(d,i);
    end
end


% yfilter  = bilateralFilter(I,I,theta(1),theta(2),theta(3),theta(4));

% 异常值修正
% h = fspecial('average',3);
% I2 = imfilter(I,h,'replicate','same');
% mask = norm01(abs(I-I2))>0.8;
% I(mask) = I2(mask);

I2 = I; % 潜在清晰图
for i =1:count
    
    
down_sample_rate4 = 1/8;
down_sample_rate3 = 1/4;
down_sample_rate2 = 1/2;
I_down4 = imresize(I,[size_I_x*down_sample_rate4, nan]);
I_down3 = imresize(I,[size_I_x*down_sample_rate3, nan]);
I_down2 = imresize(I,[size_I_x*down_sample_rate2, nan]);
% estimate params of bilateral Filter
% I,ed,step,start1,end1,start2,end2,ita,alpha,beta


range1 = [1,100];
range2 = [1,100];

step =4;
% 
[theta,d] = adaptBilateralFilter(I_down4,...
    I_down4*0,...
    step,...
    range1,...
    range2,...
    down_sample_rate4,...
    ita,...
    alpha,...
    beta);
% % 
range1 = [max(1,theta(1)-step/down_sample_rate4), min(theta(1)+step/down_sample_rate4,range1(2))];
range2 = [max(1,theta(3)-step/down_sample_rate4), min(theta(3)+step/down_sample_rate4,range2(2))];
step = 4;
% 
[theta,d] = adaptBilateralFilter(I_down3,...
    I_down3*0,...
    step,...
    range1,...
    range2,...
    down_sample_rate3,...
    ita,...
    alpha,...
    beta);


range1 = [max(1,theta(1)-step/down_sample_rate3), min(theta(1)+step/down_sample_rate3,range1(2))];
range2 = [max(1,theta(3)-step/down_sample_rate3), min(theta(3)+step/down_sample_rate3,range2(2))];

step = 4;

[theta,d] = adaptBilateralFilter(I_down2,...
    I_down2*0,...
    step,...
    range1,...
    range2,...
    down_sample_rate2,...
    ita,...
    alpha,...
    beta);

theta(1) = theta(1)/down_sample_rate2;
theta(3) = theta(3)/down_sample_rate2;
theta(2) = 2.5;


    
    yfilter = bilateralFilter(I2,I2,theta(1),theta(2),theta(3),theta(4));
    
    b= bezierFittingB(yfilter,d,C_N_I);
    %     th = mean(b(:))-mean([b(1,1),b(end,end),b(1,end),b(end,1)]);% early stop
    %     if(th < 0)
    %         break;
    %     end
    I_sub2 = I2-gamma.*b;
    
    % change weight
    gamma = max(0.01,abs(gamma-0.01)); % 权可变
    I_res2 = I_sub2 + abs(min(I_sub2(:)));
    I2 = norm01(I_res2);
    if d>2 %&& length(find(imregionalmax(yfilter))) > d
        %%%%%%% for test only
        %         change d
        d =d-1;
    end
    
    
end
I_out = I2;
end