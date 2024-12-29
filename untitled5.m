% 1. 设置参数
f_R = 40000; % 传感器谐振频率
sigma = 0.001; % 高斯脉冲的标准差，控制脉冲宽度
a = 0.02; % 接收孔径半径
c = 343; % 声速
z = 0.5; % 传播距离，假设值

% 2. 计算时间范围
t = -0.005:0.000001:0.005; % 时间范围，步长需足够小以准确表示波形
w = exp(-t.^2 / (2*sigma^2)).*sin(2*pi*f_R*t);

% 绘制发射脉冲波形w(t)
figure;
plot(t, w);
title('发射脉冲波形 w(t)');
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 设置入射角度范围
alpha = [0, 10, 20, 30, 40]; % 可选择多个角度：0° 和 20° 示例

% 初始化脉冲响应矩阵
h_R = zeros(length(t), length(alpha)); % 用于存储不同角度下的h_R
h_T_R = zeros(length(t), length(alpha)); % 用于存储不同角度下的h_T_R
r = zeros(length(t), length(alpha)); % 用于存储不同角度下的回波波形

% 设置目标幅度为 5 * 10000，并添加一个更大的幅度因子 K
target_amplitude = 5 * 10000; 
K = 1e6; % 高斯脉冲幅度因子，确保幅度远大于 5 * 10000

% 处理每个入射角度
for i = 1:length(alpha)
    angle = alpha(i) * pi / 180; % 将角度转换为弧度
    
    % 处理0°入射角（即垂直入射）
    if alpha(i) == 0
        % 计算目标时间 t = 2z/c
        target_time = 2 * z / c;

        % 找到离 target_time 最近的时间点的索引
        [~, idx] = min(abs(t - target_time));

        % 初始化 h_R 和 h_T_R 为零
        h_R(:,i) = zeros(size(t));
        h_T_R(:,i) = zeros(size(t));

        % 在目标时间位置插入一个高斯脉冲，并调整幅度
        h_R(:,i) = (K / (sigma * sqrt(2*pi))) * exp(-(t - t(idx)).^2 / (2 * sigma^2));
        h_T_R(:,i) = (K / (sigma * sqrt(2*pi))) * exp(-(t - t(idx)).^2 / (2 * sigma^2));

    else
        % 处理非零入射角
        % 计算w_squared
        w_squared = (c^2 * (t - 2*z/c).^2) / (a^2 * sin(angle)^2);

        % 限制w_squared的最大值，防止其过大
        w_squared(w_squared > 1) = 1;

        % 计算h_R（接收孔径脉冲响应）
        h_R(:,i) = (2*c*cos(angle) / (pi*a*sin(angle))) .* sqrt(1 - w_squared);

        % 限制h_R的有效范围
        h_R(t < (2*z - a*sin(angle))/c | t > (2*z + a*sin(angle))/c, i) = 0;

        % 使用h_R计算T/R对脉冲响应h_T_R
        h_T_R(:,i) = conv(w, h_R(:,i), 'same');
    end
    
    % 计算回波波形r(t)
    r(:,i) = conv(w, h_T_R(:,i), 'same'); % 使用卷积计算回波波形
    
    % 归一化回波波形 r(t)，确保幅值不超过1
end

max_r = max(abs(r), [], 'all'); % 找到所有回波的最大幅度
r = r / max_r; % 归一化所有回波波形，使得幅值不超过1

% 绘制回波波形r(t) 在同一张图中
figure;
hold on; % 开启多次绘制
for i = 1:length(alpha)
    plot(t, r(:,i)); % 绘制每个角度的回波波形
end
hold off; % 关闭多次绘制

% 添加图例
legend(arrayfun(@(x) ['入射角 = ', num2str(x), '°'], alpha, 'UniformOutput', false));
title('回波波形 r(t) 对不同入射角');
xlabel('时间 (s)');
ylabel('幅度');
grid on;
