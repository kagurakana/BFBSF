% Progressive BFBSF Algorithm 
function I_out = IR_correction(I,ita,alpha,beta,gamma)

addpath('../utils');

count =20;   % 迭代次数


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
I_down = imresize(I,[size_I_x*down_sample_rate, nan]);

% estimate params of bilateral Filter
[theta,d] = adaptBilateralFilter(I_down,I_down*0,4,ita,alpha,beta);

theta(1) = theta(1)/down_sample_rate;
theta(3) = theta(3)/down_sample_rate;
theta(2) = 2.5;

C_N_I = zeros(d);

for i = 2:d
    for j = 0:i
        C_N_I(i-1,j+1) =nchoosek(d,i);
    end
end


yfilter  = bilateralFilter(I,I,theta(1),theta(2),theta(3),theta(4));

% 异常值修正
% h = fspecial('average',3);
% I2 = imfilter(I,h,'replicate','same');
% mask = norm01(abs(I-I2))>0.8;
% I(mask) = I2(mask);

I2 = I; % 潜在清晰图
for i =1:count
    b= bezierFittingB(yfilter,d,C_N_I);
%     b(b<0)=0;
%     th = mean(b(:))-mean([b(1,1),b(end,end),b(1,end),b(end,1)]);% early stop
%     if(th < 0)
%         break;
%     end
    I_sub2 = I2-gamma.*b;
    
    % change weight
   gamma = max(0.01,abs(gamma-0.01)); % 权可变
    I_res2 = I_sub2 + abs(min(I_sub2(:)));
    I2 = norm01(I_res2);
    yfilter =bilateralFilter(I2,I2,theta(1),theta(2),theta(3),theta(4));
    

    if d>2 && length(find(imregionalmax(yfilter))) > d
%         change d
        d =d-1; 
    end
    
    
end
I_out = I2;
end