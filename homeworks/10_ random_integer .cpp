#include <iostream>
#include <random>

int main() {
    const int iterations = 10000;
    std::vector<int> counts(3, 0); // 用于统计每个值出现的次数
    srand(time(0));

    for (int i = 0; i < iterations; ++i) {
        double rand_num = (double)std::rand() / RAND_MAX; // 生成一个随机浮点数

        // 根据概率区间判断生成的整数
        int num = 0;
        if (rand_num < 0.25) {
            num = 0;
        }
        else if (rand_num < 0.75) {
            num = 1;
        }
        else {
            num = 2;
        }

        counts[num]++; // 统计每个值出现的次数
    }

    // 输出结果
    std::cout << "Counts:" << std::endl;
    std::cout << "0: " << counts[0] << std::endl;
    std::cout << "1: " << counts[1] << std::endl;
    std::cout << "2: " << counts[2] << std::endl;

    system("pause");
    return 0;
}
