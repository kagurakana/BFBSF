%
% 22/1/29
% 测试投票梯度统计约束清晰图
%
function [theta,d_out] = adaptBilateralFilter(I,ed,step,ita,alpha,beta)

if ~exist('step','var')
    step = 4;
end

if ~exist('ita','var')
    ita = 0.2;
end

if ~exist('alpha','var')
    alpha = 0.5;
end

if ~exist('beta','var')
    beta = 0.8;
end
d = 4;
C_N_I = zeros(d);

for i = 2:d
    for j = 0:i
        C_N_I(i-1,j+1) =nchoosek(d,i);
    end
end

I_center = bezierFittingB(I,d,C_N_I);

[H,D] = imhist(I_center);

[~,idx_H] = max(H);
P = imbinarize(I_center,D(idx_H)+0.1);

% --- start 5/25修改
% T = otsuthresh(H);
% [~,idx_H] = max(H);
% P = imbinarize(I_center,T);
% --- end 5/25修改
[I_center_p,center_idx_i_Mat] = max(I_center);
[~,center_idx_j] = max(I_center_p);

center_idx_i = center_idx_i_Mat(center_idx_j);



bias_percentage = sum(P(:))/(size(I,1)*size(I,2));
if 0.2 < bias_percentage && bias_percentage < 0.3
    range = 13;
    i_count = 13;
    j_count = 13;
%     alpha = 0.7;
    d_out = 10;
elseif sum(P(:))/(size(I,1)*size(I,2)) < 0.2
    range = 6;
    i_count = 6;
    j_count = 6;
%     alpha = 0.7;
    d_out = 15;
else 
    range = 10;
    i_count = 20;
    j_count = 20;
    d_out =6; %7->6
end

R_B = zeros(i_count);

z_f_b = zeros(i_count);
R_F = zeros(i_count);


I_iter = I;

for i=1:i_count
    for j =1:j_count
        b =bilateralFilter(I_iter,ed,(i)*step,2,(j)*step,0.7);
        
%         b = b-min(b(:));
%         [i,j]
        I_sub2 = I_iter-b;
        I2_stretching = norm01(I_sub2);
        
        if sum(P(:))~=0
            I2_mask = I2_stretching.*P;
        else
            I2_mask = I2_stretching;
        end
        
        I_00 = I2_mask(1:center_idx_i,1:center_idx_j);
        I_01 = I2_mask(1:center_idx_i,center_idx_j+1:end);
        I_10 = I2_mask(center_idx_i+1:end,1:center_idx_j);
        I_11 = I2_mask(center_idx_i+1:end,center_idx_j+1:end);
        
        if sum(P(:))~=0
            b_mask = b.*P;
        else
            b_mask = b;
        end
        
        b_00 = b_mask(1:center_idx_i,1:center_idx_j);
        b_01 = b_mask(1:center_idx_i,center_idx_j+1:end);
        b_10 = b_mask(center_idx_i+1:end,1:center_idx_j);
        b_11 = b_mask(center_idx_i+1:end,center_idx_j+1:end);
        
        % countHOG(Image, filter, bar_count, is_weight)
        %%%%%
        %   ↖ ↑↗     2  |  1
        %   ← O →    --  9 --
        %   ↙↓ ↘     4  |  3
        %%%%
       
        b_00_HOG = countHOG(b_00,6,4,1);
        b_01_HOG = countHOG(b_01,6,4,1);
        b_10_HOG = countHOG(b_10,6,4,1);
        b_11_HOG = countHOG(b_11,6,4,1);
        D_B = imgradient(b,'roberts');
        R_B(i,j) = sum(abs(D_B(:)))/(b_00_HOG(3)+...
            b_01_HOG(4)+...
            b_10_HOG(1)+...
            b_11_HOG(2));            % BTV
        % 块方向
        
        I_00_HOG = voteCountHOG(I_00,6,4);
        I_01_HOG = voteCountHOG(I_01,6,4);
        I_10_HOG = voteCountHOG(I_10,6,4);
        I_11_HOG = voteCountHOG(I_11,6,4);
        
        mean_I_00 = mean(I_00_HOG);
        mean_I_01 = mean(I_01_HOG);
        mean_I_10 = mean(I_10_HOG);
        mean_I_11 = mean(I_11_HOG);
        
        
        R_F(i,j) = sum([...
            norm(I_00_HOG-mean_I_00,2),...
            norm(I_01_HOG-mean_I_01,2),...
            norm(I_10_HOG-mean_I_10,2),...
            norm(I_11_HOG-mean_I_11,2),...
            ]);
            
        z_f_b(i,j) = norm(reshape(I-I2_stretching-b,1,[]),2);
        
%         t=[i,j]
    end
end

R_B = featureScaling(R_B);
R_F = featureScaling(R_F);
z_f_b = featureScaling(z_f_b);
% model
constraints =  ita.*z_f_b + alpha.*R_F + beta.*R_B;
constraints = constraints(1:range,1:range);
% -Rb
% constraints =  ita.*z_f_b+alpha.*mean_dir_I;
% -Rf
% constraints =  ita.*z_f_b+beta.*db_grad;

[min_d,idx_ii] = min(constraints);
[~,idx_jj] = min(min_d);

ii = idx_ii(idx_jj);
jj = idx_jj;

% figure
% mesh(constraints);

theta = [(ii)*step,2,(jj)*step,0.7];

end
function mormlizationMat = featureScaling(Mat)
mormlizationMat = (Mat-min(Mat(:)))./(max(Mat(:)) - min(Mat(:)))+1;
% mormlizationMat = Mat./max(Mat(:));
% mormlizationMat = mormlizationMat-mean(mormlizationMat(:));
end
