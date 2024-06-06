#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <corecrt_math_defines.h>


const double g = 9.81; // 重力加速度

// 双摆系统参数
struct Pendulum {
    double m1, m2; // 质量
    double l1, l2; // 长度
};

// 定义状态向量
struct State {
    double theta1, omega1, theta2, omega2;
};

// 双摆系统的微分方程
std::vector<double> equations(double t, const State& state, const Pendulum& pendulum) {
    double delta_theta = state.theta2 - state.theta1;

    double den1 = (pendulum.m1 + pendulum.m2) * pendulum.l1 - pendulum.m2 * pendulum.l1 * cos(delta_theta) * cos(delta_theta);
    double den2 = (pendulum.l2 / pendulum.l1) * den1;

    double dtheta1_dt = state.omega1;
    double dtheta2_dt = state.omega2;
    double domega1_dt = (pendulum.m2 * pendulum.l1 * state.omega1 * state.omega1 * sin(delta_theta) * cos(delta_theta) +
        pendulum.m2 * g * sin(state.theta2) * cos(delta_theta) +
        pendulum.m2 * pendulum.l2 * state.omega2 * state.omega2 * sin(delta_theta) -
        (pendulum.m1 + pendulum.m2) * g * sin(state.theta1)) / den1;
    double domega2_dt = (-pendulum.m2 * pendulum.l2 * state.omega2 * state.omega2 * sin(delta_theta) * cos(delta_theta) +
        (pendulum.m1 + pendulum.m2) * (g * sin(state.theta1) * cos(delta_theta) -
            pendulum.l1 * state.omega1 * state.omega1 * sin(delta_theta) - g * sin(state.theta2))) / den2;

    return { dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt };
}

// 四阶Runge-Kutta方法求解器
State runge_kutta_step(double t, const State& state, double dt, const Pendulum& pendulum) {
    auto k1 = equations(t, state, pendulum);
    State state_k2 = { state.theta1 + 0.5 * dt * k1[0], state.omega1 + 0.5 * dt * k1[1],
                       state.theta2 + 0.5 * dt * k1[2], state.omega2 + 0.5 * dt * k1[3] };
    auto k2 = equations(t + 0.5 * dt, state_k2, pendulum);

    State state_k3 = { state.theta1 + 0.5 * dt * k2[0], state.omega1 + 0.5 * dt * k2[1],
                       state.theta2 + 0.5 * dt * k2[2], state.omega2 + 0.5 * dt * k2[3] };
    auto k3 = equations(t + 0.5 * dt, state_k3, pendulum);

    State state_k4 = { state.theta1 + dt * k3[0], state.omega1 + dt * k3[1],
                       state.theta2 + dt * k3[2], state.omega2 + dt * k3[3] };
    auto k4 = equations(t + dt, state_k4, pendulum);

    State new_state;
    new_state.theta1 = state.theta1 + (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    new_state.omega1 = state.omega1 + (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    new_state.theta2 = state.theta2 + (dt / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    new_state.omega2 = state.omega2 + (dt / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);

    return new_state;
}

double calculate_lyapunov_exponent(double delta_0, const std::vector<double>& times, const std::vector<double>& theta1_1, const std::vector<double>& theta1_2) {
    std::vector<double> log_deltas(times.size());
    for (size_t i = 0; i < times.size(); ++i) {
        double delta_t = std::abs(theta1_1[i] - theta1_2[i]);
        log_deltas[i] = std::log(delta_t / delta_0);
    }

    // 使用线性拟合计算Lyapunov指数
    double sum_t = 0.0, sum_log_delta = 0.0, sum_t_log_delta = 0.0, sum_t2 = 0.0;
    for (size_t i = 0; i < times.size(); ++i) {
        sum_t += times[i];
        sum_log_delta += log_deltas[i];
        sum_t_log_delta += times[i] * log_deltas[i];
        sum_t2 += times[i] * times[i];
    }

    double n = times.size();
    double lyapunov_exponent = (n * sum_t_log_delta - sum_t * sum_log_delta) / (n * sum_t2 - sum_t * sum_t);
    return lyapunov_exponent;
}


int main() {
    // 初始条件////////////修改初始条件在此处/////////////
    State initial_state = { M_PI / 2, 0.0, M_PI / 2, 0.0 };
    State perturbed_state = { M_PI / 2 + 1e-8, 0.0, M_PI / 2, 0.0 };


    // 系统参数
    Pendulum pendulum = { 1.0, 1.0, 1.0, 1.0 };

    // 时间参数
    double t_start = 0.0;
    double t_end = 10.0;
    double dt = 0.01;
    size_t num_steps = static_cast<size_t>((t_end - t_start) / dt);

    // 内存保存时间和状态数据
    std::vector<double> times(num_steps);
    std::vector<double> theta1_1(num_steps), theta1_2(num_steps);
    // 文件保存数据
    std::ofstream file("double_pendulum_data.csv");
    file << "time,theta1,omega1,theta2,omega2\n";

    // 循环求解
    State state = initial_state;
    State perturbed = perturbed_state;
    /*for (double t = t_start; t <= t_end; t += dt) {
        file << t << "," << state.theta1 << "," << state.omega1 << "," << state.theta2 << "," << state.omega2 << "\n";
        state = runge_kutta_step(t, state, dt, pendulum);
    }*/
    for (size_t i = 0; i < num_steps; ++i) {
        double t = t_start + i * dt;
        times[i] = t;
        theta1_1[i] = state.theta1;
        theta1_2[i] = perturbed.theta1;
        file << t << "," << state.theta1 << "," << state.omega1 << "," << state.theta2 << "," << state.omega2 << "\n";
        state = runge_kutta_step(t, state, dt, pendulum);
        perturbed = runge_kutta_step(t, perturbed, dt, pendulum);
    }
   
    
    // 计算Lyapunov指数
        double delta_0 = 1e-8;
    double lyapunov_exponent = calculate_lyapunov_exponent(delta_0, times, theta1_1, theta1_2);
    std::cout << "Lyapunov Exponent: " << lyapunov_exponent << std::endl;


    file.close();
    std::cout << "Simulation completed. Data saved to double_pendulum_data.csv" << std::endl;
    return 0;
}
