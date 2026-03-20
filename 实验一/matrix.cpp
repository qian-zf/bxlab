#include <iostream>
#include <windows.h>  // Sleep, GetTickCount, timeGetTime, QueryPerformanceCounter

using namespace std;

void init(int n, double* sum, double* a, double** b) {
    for (int i = 0; i < n; i++) {
        sum[i] = 0;
        a[i] = i;  // 修正 i++ 问题
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i][j] = i + j;
        }
    }
}

void pingfan(int n, double* sum, double* a, double** b) {
    for (int i = 0; i < n; i++) {
        sum[i] = 0;
        for (int j = 0; j < n; j++)
            sum[i] += b[j][i] * a[j];
    }
}

void youhua(int n, double* sum, double* a, double** b) {
    for (int i = 0; i < n; i++)
        sum[i] = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            sum[i] += b[j][i] * a[j];
        }
    }
}

// 计时器封装函数
double benchmark(void (*func)(int, double*, double*, double**), int n, double* sum, double* a, double** b, int counter) {
    double total_time = 0;
    for (int i = 0; i < counter; i++) {
        LARGE_INTEGER litmp;
        LONGLONG QPart1, QPart2;
        double dfMinus, dfFreq, dfTim;
        QueryPerformanceFrequency(&litmp);
        dfFreq = (double)litmp.QuadPart;  // 获得计数器的时钟频率
        QueryPerformanceCounter(&litmp);
        QPart1 = litmp.QuadPart;  // 获得初始值

        func(n, sum, a, b);  // 调用传入的算法函数

        QueryPerformanceCounter(&litmp);
        QPart2 = litmp.QuadPart;  // 获得终止值
        dfMinus = (double)(QPart2 - QPart1);
        dfTim = dfMinus * 1000 / dfFreq;
        total_time += dfTim;
    }
    return total_time / counter;
}

int main() {
    // 批量测试不同矩阵大小
    for (int n = 100; n <= 4000; n += 100) {
        double* sum = new double[n];
        double* a = new double[n];
        double** b = new double*[n];
        for (int i = 0; i < n; i++) {
            b[i] = new double[n];
        }

        init(n, sum, a, b);

        int counter = 40;
        double avg_time = benchmark(pingfan, n, sum, a, b, counter);
        cout << "n = " << n << ", 平均时间: " << avg_time << " ms" << endl;

        // 释放内存
        delete[] sum;
        delete[] a;
        for (int i = 0; i < n; i++) {
            delete[] b[i];
        }
        delete[] b;
    }

    // 手动输入矩阵大小和测试次数
    int n, counter;
    cout << "输入矩阵大小: ";
    cin >> n;
    cout << "输入测试次数: ";
    cin >> counter;

    double* sum = new double[n];
    double* a = new double[n];
    double** b = new double*[n];
    for (int i = 0; i < n; i++) {
        b[i] = new double[n];
    }
    init(n, sum, a, b);

    double avg_time = benchmark(youhua, n, sum, a, b, counter);
    cout << "优化算法平均时间: " << avg_time << " ms" << endl;

    // 释放内存
    delete[] sum;
    delete[] a;
    for (int i = 0; i < n; i++) {
        delete[] b[i];
    }
    delete[] b;

    return 0;
}
