function bias = Bias_gen_norm(img,brightness,radii,center_x,center_y,min_val,max_val)
% 生成热辐射退化图

img = im2gray(img);
img = im2double(img);
[Width,Height,dim]=size(img);
mu_w=Width*center_x;
mu_h=Height*center_y;                   %位置
sigma =min(Width,Height)*radii;      %方差 控制圈大小

%生成高斯曲面
u = 0:1:(Height -1);
v = 0:1:(Width -1);
[U,V] = meshgrid(u,v);
H=exp(-((U-mu_h).^2+(V-mu_w).^2)./(2*sigma^2))*brightness;% 亮度
bias = norm_min_max(H,min_val,max_val);
end