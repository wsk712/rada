% 参数设置
f_R = 40000;          % 传感器谐振频率 (单位: Hz)
sigma = 0.001;        % 脉冲持续时间相关参数 (单位: s)
a = 0.02;             % 接收孔径半径 (单位: m)
c = 343;              % 声速 (单位: m/s)
z = 0.5;              % 传播距离 (单位: m)
theta_deg = 20;       % 入射角 (单位: 度)
theta = theta_deg * pi / 180; % 转换为弧度

% 计算时间范围
t = -0.01:0.000001:0.01; % 时间范围，步长为 1 µs
dt = t(2) - t(1);        % 时间采样间隔，用于卷积后归一化
w = exp(-t.^2 / (2 * sigma^2)) .* sin(2 * pi * f_R * t); % 发射脉冲波形

% 绘制发射脉冲波形 w(t)
figure;
plot(t, w);
title('发射脉冲波形 w(t)');
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 计算接收孔径脉冲响应 h_R(t)
h_R = zeros(size(t)); % 初始化 h_R
w_squared = (c^2 * (t - 2 * z / c).^2) / (a^2 * sin(theta)^2); % w_squared 计算

% 确定有效时间范围
valid_time_range = (t >= (2 * z - a * sin(theta)) / c) & ...
                   (t <= (2 * z + a * sin(theta)) / c);

% 在有效时间范围内计算 h_R(t)
h_R(valid_time_range) = (2 * c * cos(theta) / (pi * a * sin(theta))) .* ...
                        sqrt(1 - w_squared(valid_time_range));

% 绘制 h_R(t) 波形
figure;
plot(t, h_R);
title(['接收孔径脉冲响应 h_R(t) (入射角 = ', num2str(theta_deg), '°)']);
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 使用卷积计算 h_T_R(t)
h_T_R = conv(h_R, h_R, 'same') * dt; % 卷积结果并乘以采样间隔

% 绘制 h_T_R(t) 波形
figure;
plot(t, h_T_R);
title(['接收孔径脉冲响应 h_{T/R}(t) (使用卷积计算，入射角 = ', num2str(theta_deg), '°)']);
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 计算最终响应 r(t) = h_T_R(t) * w(t) (卷积)
r_t = conv(h_T_R, w, 'same') * dt; % 卷积结果并乘以采样间隔

% 绘制 r(t) 波形
figure;
plot(t, r_t);
title('最终响应 r(t) = h_{T/R}(t) \ast w(t)');
xlabel('时间 (s)');
ylabel('幅度');
grid on;
