% 参数设置
f_R = 40000; % 传感器谐振频率
sigma = 0.001; % 脉冲持续时间相关参数
a = 0.02; % 接收孔径半径
c = 343; % 声速
z = 0.5; % 传播距离，假设值

% 计算时间范围
t = -0.01:0.000001:0.01; % 时间范围，步长需足够小以准确表示波形
w = exp(-t.^2 / (2*sigma^2)).*sin(2*pi*f_R*t);

% 绘制发射脉冲波形w(t)
figure;
plot(t, w);
title('发射脉冲波形 w(t)');
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 入射角20°
theta_deg = 40; % 角度，单位是度
theta = theta_deg * pi / 180; % 转换为弧度

% 计算w_squared
w_squared = (c^2 * (t - 2*z/c).^2) / (a^2 * sin(theta)^2);

% 计算h_R
h_R = zeros(size(t)); % 初始化h_R为零

% 计算有效时间范围
valid_time_range = (t >= (2*z - a*sin(theta))/c) & (t <= (2*z + a*sin(theta))/c);

% 仅在有效时间范围内计算h_R
h_R(valid_time_range) = (2*c*cos(theta) / (pi*a*sin(theta))) .* sqrt(1 - w_squared(valid_time_range));

% 绘制h_R波形
figure;
plot(t, h_R);
title(['接收孔径脉冲响应 h_R(t) (入射角 = ', num2str(theta_deg), '°)']);
xlabel('时间 (s)');
ylabel('幅度');
grid on;

% 定义 h_T_R 的计算
h_T_R = zeros(size(t)); % 初始化 h_T_R

% 使用数值积分计算 h_T_R
for i = 1:length(t)
    % 当前时间点 t(i)
    t_current = t(i);

    % 定义当前时间点下的被积函数
    integrand = @(tau) interp1(t, h_R, tau, 'linear', 0) .* ...
                       interp1(t, h_R, t_current - tau, 'linear', 0);

    % 积分上下限
    lower_limit = (2*z - a*sin(theta)) / c;
    upper_limit = (2*z + a*sin(theta)) / c;

    % 检查积分上下限是否在范围内
    if lower_limit >= min(t) && upper_limit <= max(t)
        % 计算数值积分
        h_T_R(i) = integral(integrand, lower_limit, upper_limit, 'ArrayValued', true);
    else
        h_T_R(i) = 0; % 超出范围，设为0
    end
end

% 绘制 h_T_R 波形
figure;
plot(t, h_T_R);
title(['接收孔径脉冲响应 h_{T/R}(t) (入射角 = ', num2str(theta_deg), '°)']);
xlabel('时间 (s)');
ylabel('幅度');
grid on;
