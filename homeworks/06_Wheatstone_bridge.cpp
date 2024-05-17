#define _CRT_SECURE_NO_WARNINGS //���԰�ȫ���
#include <iostream>
#include <cmath>

#define PI 3.14159265358

// ��ʼ������
void initializeMatrix(double matrix[3][3], const double r1, const double r2, const double r3, const double rx, const double ra, const double rs)
{
    matrix[0][0] = rs;
    matrix[0][1] = r1;
    matrix[0][2] = r2;
    matrix[1][0] = -rx;
    matrix[1][1] = r1 + rx + ra;
    matrix[1][2] = -ra;
    matrix[2][0] = -r3;
    matrix[2][1] = -ra;
    matrix[2][2] = r2 + r3 + ra;
}

void showMatrix(double matrix[3][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl << std::endl << std::endl;
    }
}
//����ʽ
double Determinant(double matrix[3][3])
{
    double det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    return det;
}
//�������
void Adjoint(double matrix[3][3], double adj[3][3]) 
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            // �����Ӿ�������ʽ
            double subDet = matrix[(i + 1) % 3][(j + 1) % 3] * matrix[(i + 2) % 3][(j + 2) % 3] -
                matrix[(i + 1) % 3][(j + 2) % 3] * matrix[(i + 2) % 3][(j + 1) % 3];
            // ����λ�ü����������Ԫ��
            adj[j][i] = subDet;
        }
    }
}
//�����
bool Inverse(double matrix[3][3], double detValue, double inverse[3][3]) 
{
    // �������ʽΪ0��������󲻴���
    if (detValue == 0) {
        return false;
    }

    double adjoint[3][3];
    Adjoint(matrix, adjoint);

    // ���������
    double invDet = 1.0 / detValue;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            inverse[i][j] = adjoint[i][j] * invDet;
        }
    }

    return true;
}
//[3*3]X[3*1]
void Multiply(double matrix[3][3], double vector[3], double result[3]) {
    for (int i = 0; i < 3; ++i) {
        result[i] = 0;
        for (int j = 0; j < 3; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

int main()
{

	const int MATRIX_SIZE = 3;
    double r1, r2, r3, rx, ra, rs;
    double i[3] = { 0 ,0 ,0 };
    double v[3] = { 0 ,0 ,0 };
	std::cout << "��������� r1, r2, r3, rx, ra, rs �͵�Դ��ѹ v0." << std::endl;
    std::cin >> r1 >> r2 >> r3 >> rx >> ra >> rs >> v[0];

    double Resistance[3][3];
    initializeMatrix(Resistance, r1, r2, r3, rx, ra, rs);
    std::cout << "����Ϊ��" << std::endl;
    showMatrix(Resistance);

    double determinant = Determinant(Resistance);
    double inverseMatrix[3][3];

    bool success = Inverse(Resistance, determinant, inverseMatrix);

    if (success) {
        std::cout << "�����Ϊ��" << std::endl;
        showMatrix(inverseMatrix);
    }
    else {
        std::cout << "���������ʽΪ0������󲻴��ڣ��޷���⣡" << std::endl;
        return 1;
    }
    Multiply(inverseMatrix, v, i);
    std::cout << "The equivalent resistance:" << v[0] / i[0] << std::endl;
    system("pause");
    return 0;
}