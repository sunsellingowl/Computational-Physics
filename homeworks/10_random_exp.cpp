#include <iostream>
#include <random>

// Define the distribution function
double distribution(double x) {
    return exp(-pow(x - 1, 2));
}

double generateRandomNumber() {
    std::random_device rd;
    std::mt19937 gen(rd());

    double x, y;
    do {
        // Generate random x and y values
        std::uniform_real_distribution<> dis_x(0, 2);
        std::uniform_real_distribution<> dis_y(0, 1);
        x = dis_x(gen);
        y = dis_y(gen);

        // Check if y is less than the probability density function
    } while (y > distribution(x));

    return x;
}

int main() {
    const int iterations = 100;
    // Generate and print 100 random numbers
    std::cout << "Generated random numbers with the given distribution:\n";
    for (int i = 0; i < iterations; ++i) {
        std::cout << generateRandomNumber() << '\t';
        if ((i + 1) % 5 == 0) {
            std::cout << std::endl;
        }
    }
    system("pause");
    return 0;
}