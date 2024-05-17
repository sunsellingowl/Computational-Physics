#include <iostream>
#include <random>

int main() {
    const int iterations = 10000;
    std::vector<int> counts(3, 0); // ����ͳ��ÿ��ֵ���ֵĴ���
    srand(time(0));

    for (int i = 0; i < iterations; ++i) {
        double rand_num = (double)std::rand() / RAND_MAX; // ����һ�����������

        // ���ݸ��������ж����ɵ�����
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

        counts[num]++; // ͳ��ÿ��ֵ���ֵĴ���
    }

    // ������
    std::cout << "Counts:" << std::endl;
    std::cout << "0: " << counts[0] << std::endl;
    std::cout << "1: " << counts[1] << std::endl;
    std::cout << "2: " << counts[2] << std::endl;

    system("pause");
    return 0;
}
