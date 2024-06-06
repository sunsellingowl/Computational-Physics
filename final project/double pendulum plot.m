data = readtable(['E:\Projects\MSVS\Computational Physics\double pendulum' ...
    '\double pendulum\double_pendulum_data.csv']);

% 提取时间和角度数据
time = data.time;
theta1 = data.theta1;
theta2 = data.theta2;

% 双摆的长度
l1 = 1.0;
l2 = 1.0;

% 绘制角度随时间的变化
% figure;
% plot(time, theta1, 'r', 'DisplayName', '\theta_1');
% hold on;
% plot(time, theta2, 'b', 'DisplayName', '\theta_2');
% xlabel('Time [s]');
% ylabel('Angle [rad]');
% legend;
% title('Double Pendulum Angles');
% grid on;

% 计算双摆末端的坐标
x1 = l1 * sin(theta1);
y1 = -l1 * cos(theta1);
x2 = x1 + l2 * sin(theta2);
y2 = y1 - l2 * cos(theta2);

initconditions = "\theta_1 = \pi/2, \theta_2 = -\pi/2";
% 绘制双摆的运动轨迹
figure;
plot(x1, y1, 'r', 'DisplayName', 'Pendulum 1','LineWidth',1.1);
hold on;
plot(x2, y2, 'b', 'DisplayName', 'Pendulum 2','LineWidth',1.1);
xlabel('X');
ylabel('Y');
legend;
title('Double Pendulum Trajectory');
text(-0.4,0.4,initconditions);
grid on;



