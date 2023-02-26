%
function b = bezierFittingB(yfilter,d,C_N_I)
    % b: 估计热辐射曲面
    % I: 退化图
    % d: 曲面次数  
%     [xx,yy]=meshgrid(1:size(yfilter,2),1:size(yfilter,1));
%     meshc(xx,yy,yfilter);title('双边滤波图像');
    %%
 
    row = ceil(size(yfilter,2)/32);
    col = ceil(size(yfilter,1)/32);

% --- start 5/25 修改
%     row = ceil(size(yfilter,2)/128);
%     col = ceil(size(yfilter,1)/128);
   
% --- end 5/25 修改
%     window_size = [32,32];
 window_size = [row,col];
    % 下采样
    I2 = blockproc(yfilter, window_size, @block_mid);
%     imshow(I2);

    [given_point_x,given_point_y] = meshgrid(...
        1:window_size(2):(size(yfilter,2)),...
        1:window_size(1):(size(yfilter,1))...
        );

    given_point_z = I2;


    %%
    given_point_count = size(given_point_x,1) * size(given_point_x,2);
    if ~exist('d','var')
        m = 7; % u_deriction
        n = 7; % v_deriction
    else
        m = d;
        n = d;
    end
    flat_given_x = reshape(given_point_x,[1,given_point_count]);
    flat_given_y = reshape(given_point_y,[1,given_point_count]);
    flat_given_z = reshape(given_point_z,[1,given_point_count]);

    % u 方向求差分 ↓
    diff_u = sqrt(diff(given_point_x,1,1).^2+...
        diff(given_point_y,1,1).^2+...
        diff(given_point_z,1,1).^2);
    % v 方向求差分 →
    diff_v = sqrt(diff(given_point_x,1,2).^2+...
        diff(given_point_y,1,2).^2+...
        diff(given_point_z,1,2).^2);

    percentage_p_u = diff_u ./ sum(diff_u);
    percentage_p_v = diff_v ./ sum(diff_v,2);

    u = zeros(size(given_point_x,1), size(given_point_x,2));
    v = zeros(size(given_point_x,1), size(given_point_x,2));

    for i =2:size(u,1)
        u(i, :) = u(i-1, :) + percentage_p_u(i-1, :);
    end

    for i =2:size(v,2)
        v(:,i) = v(:, i-1) + percentage_p_v(:, i-1);
    end
  
    % 有多少点就会有多少行
    Bernstein_mat = zeros(given_point_count,(m+1)*(n+1));

    flat_u = reshape(u,[1,given_point_count]);
    flat_v = reshape(v,[1,given_point_count]);
    for i =1:size(Bernstein_mat,1)
        for j = 0:m % B(u)
            temp = B_n_i(m, j, flat_u(i),C_N_I);
            for k = 0:n %B(v)
                % B_n_i(m, j, flat_u(i))*B_n_i(n, k,flat_v(i));
                Bernstein_mat(i,j*(n+1)+k+1) = temp*B_n_i(n, k,flat_v(i),C_N_I);
            end
        end
    end


     control_x = lsqlin(Bernstein_mat,flat_given_x');
     control_y = lsqlin(Bernstein_mat,flat_given_y');
     control_z = lsqlin(Bernstein_mat,flat_given_z');



     % ======= PLOT =========

     control_x = reshape(control_x,[m+1,n+1]);
     control_y = reshape(control_y,[m+1,n+1]);
    control_z = reshape(control_z,[m+1,n+1]);

    size_x = size(control_z,1); %x方向控制点个数
    size_y = size(control_z,2); %y方向控制点个数



%      flat_control_x = reshape(control_x,[1,control_point_count]);
%      flat_control_y = reshape(control_y,[1,control_point_count]);
%      flat_control_z = reshape(control_z,[1,control_point_count]);


    %point_count_x = 40; % 精度，两个方向的取点数
    %point_count_y = 40;

    point_count_x = size(yfilter,1); % 精度，两个方向的取点数
    point_count_y = size(yfilter,2);


    p_x = zeros(point_count_x,point_count_y);p_y = zeros(point_count_x,point_count_y);p_z = zeros(point_count_x,point_count_y);

    [t_u, t_v] = meshgrid(linspace(0,1,point_count_y), linspace(0,1,point_count_x)); % 生成u_v平面




    for i = 1:size_x
        bernstein_u = B_n_i(size_x-1, i-1, t_u,C_N_I);% 第 i - 1 个Bernstein基
        for j = 1:size_y
            bernstein_v = B_n_i(size_y-1, j-1, t_v,C_N_I);% 第 j - 1 个Bernstein基
%             p_x = p_x + control_x(i,j)*bernstein_u.*bernstein_v;
%             p_y = p_y + control_y(i,j)*bernstein_u.*bernstein_v;
            p_z = p_z + control_z(i,j)*bernstein_u.*bernstein_v;
        end
    end

 
    b = p_z;

    % Bernstein 基
    function B_n_i = B_n_i(n,i,t,C_N_I)
        B_n_i = C_N_I(n-1,i+1) * t.^ (i) .* (1-t) .^ (n - i); % 注意用点乘！
    end

    % 邻域求平均
    function x = block_mid(struct)
        idx1 = ceil(struct.blockSize(1)/2);
        idx2 = ceil(struct.blockSize(2)/2);
        x = struct.data(idx1,idx2);
    end


end