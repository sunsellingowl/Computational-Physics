## 统计物理中的数值计算初步 - 详细学习笔记

### 目录
1. [Monte Carlo 初步](#monte-carlo-初步)
2. [Metropolis 算法](#metropolis-算法)
3. [2D Ising 模型的 Monte Carlo 方法](#2d-ising-模型的-monte-carlo-方法)
4. [高级算法](#高级算法)
5. [Quantum Monte Carlo 简介](#quantum-monte-carlo-简介)
6. [总结](#总结)
7. [参考文献](#参考文献)

### Monte Carlo 初步
#### 什么是 Monte Carlo 方法？
- Monte Carlo 方法是一种通过使用随机数来解决复杂数学和物理问题的数值计算方法。
- 最早在二战期间的“曼哈顿计划”中被用来模拟中子扩散，由数学家冯·诺伊曼命名。
- 主要思想是利用大量随机数进行实验或计算，以获得近似解。

#### Monte Carlo 方法的基本步骤
1. **生成随机样本**：在给定的域中生成大量随机样本点。
2. **计算函数值**：对于每一个随机样本点，计算目标函数值。
3. **求平均值**：对所有样本点的函数值取平均，作为积分或期望值的近似。
4. **估计误差**：通过样本的标准偏差估计计算结果的误差。

#### Monte Carlo 方法的特点
- **易于处理复杂问题**：适用于高维空间和复杂边界的积分问题。
- **简单的误差估计**：通过计算样本的标准偏差可以直接估计误差。
- **适用于多种应用**：广泛应用于物理、金融、工程等领域的数值模拟和优化问题。
- **收敛速度慢**：与样本数量的平方根成反比，收敛速度较慢。

#### 示例：计算 π 的值
1. 在边长为1的正方形内随机生成 \(N\) 个点。
2. 计算这些点中落在半径为1的四分之一圆内的点数 \(M\)。
3. 使用比例 \(\frac{M}{N} \approx \frac{\pi}{4}\) 估算 π 的值。


```python
import random

N = 1000000
M = 0

for _ in range(N):
    x, y = random.uniform(0, 1), random.uniform(0, 1)
    if x**2 + y**2 <= 1:
        M += 1

pi_estimate = 4 * M / N
print(f"Estimated π value: {pi_estimate}")
```


### Metropolis 算法
#### 经典统计物理中的应用
- **哈密顿量 \( H \)**：系统的能量函数，描述系统的状态。
- **配分函数 \( Z \)**：用于计算系统的宏观物理量，如自由能 \( F = -k_B T \ln Z \)。
- 直接计算配分函数涉及高维积分，计算复杂度高，适合用 Monte Carlo 方法解决。

#### Metropolis 算法步骤
1. **初始化**：选择系统的初始状态 \( x_0 \)。
2. **生成候选状态**：从当前状态 \( x_t \) 生成一个新状态 \( x_{t+1} \)。
3. **计算能量差**：计算新旧状态的能量差 \( \Delta E = E(x_{t+1}) - E(x_t) \)。
4. **接受准则**：
    - 若 \( \Delta E \leq 0 \)，接受新状态。
    - 若 \( \Delta E > 0 \)，以概率 \( e^{-\Delta E / k_B T} \) 接受新状态。
5. **重复**：不断重复以上步骤，直到系统达到平衡态。

#### 伪代码

```python
import random
import math

def metropolis_algorithm(energy, state, temperature):
    current_state = state
    current_energy = energy(current_state)
    
    while True:
        new_state = generate_new_state(current_state)
        new_energy = energy(new_state)
        delta_E = new_energy - current_energy
        
        if delta_E <= 0 or random.uniform(0, 1) < math.exp(-delta_E / temperature):
            current_state = new_state
            current_energy = new_energy
        
        yield current_state, current_energy
```


### 2D Ising 模型的 Monte Carlo 方法
- **Ising 模型**：描述自旋系统的简单模型，用于研究磁性相变。
- **2D Ising 模型**：考虑二维正方格子中的自旋，取值 \( \sigma_i = \pm 1 \)。
- **哈密顿量**：\( H = -J \sum_{\langle i,j \rangle} \sigma_i \sigma_j \)，其中 \( J \) 为交换常数， \( \langle i,j \rangle \) 表示最近邻自旋对。

#### Metropolis 算法应用于 2D Ising 模型
1. **初始化**：设置所有自旋的初始状态（例如，随机分布）。
2. **更新自旋**：随机选择一个自旋，计算翻转该自旋前后的能量差 \( \Delta E \)。
3. **接受准则**：根据 Metropolis 准则决定是否接受自旋翻转。
4. **迭代**：重复上述过程，直到系统达到平衡。

#### Python 实现示例

```python
import numpy as np
import random

def initialize_lattice(N):
    return np.random.choice([-1, 1], size=(N, N))

def calculate_energy(lattice, J=1):
    energy = 0
    N = lattice.shape[0]
    for i in range(N):
        for j in range(N):
            S = lattice[i,j]
            neighbors = lattice[(i+1)%N,j] + lattice[i,(j+1)%N] + lattice[(i-1)%N,j] + lattice[i,(j-1)%N]
            energy += -J * S * neighbors
    return energy / 2

def metropolis_step(lattice, T):
    N = lattice.shape[0]
    for _ in range(N**2):
        i, j = random.randint(0, N-1), random.randint(0, N-1)
        S = lattice[i, j]
        neighbors = lattice[(i+1)%N, j] + lattice[i, (j+1)%N] + lattice[(i-1)%N, j] + lattice[i, (j-1)%N]
        delta_E = 2 * S * neighbors
        if delta_E < 0 or random.random() < np.exp(-delta_E / T):
            lattice[i, j] *= -1

N = 20
lattice = initialize_lattice(N)
T = 2.0

for _ in range(1000):
    metropolis_step(lattice, T)

print(f"Final energy: {calculate_energy(lattice)}")
```


### 高级算法
- **模拟退火**：通过逐步降低温度，找到系统的最优解。
- **自适应采样**：根据样本分布动态调整采样策略，提高效率。
- **Cluster 算法**：如 Swendsen-Wang 算法，用于减少自旋系统的自相关时间。

### Quantum Monte Carlo 简介
- **Quantum Monte Carlo**：用于模拟量子系统的 Monte Carlo 方法，通过随机数模拟量子态的演化。
- **路径积分 Monte Carlo**：利用 Feynman 路径积分形式，计算量子系统的物理量。
- **变分 Monte Carlo**：利用变分原理，通过 Monte Carlo 积分求解基态波函数的期望值。

#### 示例：路径积分 Monte Carlo

```python
import numpy as np

def path_integral_monte_carlo(beta, N, potential):
    dtau = beta / N
    paths = np.zeros((N,))
    energies = np.zeros((N,))
    
    for i in range(1, N):
        paths[i] = paths[i-1] + np.random.normal(0, np.sqrt(dtau))
    
    for i in range(N):
        kinetic = 0.5 * (paths[(i+1)%N] - paths[i])**2 / dtau
        potential_energy = potential(paths[i])
        energies[i] = kinetic + potential_energy
    
    return np.mean(energies)

def harmonic_potential(x):
    return 0.5 * x**2

beta = 10.0
N = 100
average_energy = path_integral_monte_carlo(beta, N, harmonic_potential)
print(f"Average energy: {average_energy}")
```


