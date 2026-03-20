//Windows
#include<iostream>
#include<time.h>
#include<thread>
#include <windows.h>    //Sleep, GetTickCount, timeGetTime, QueryPerformanceCounter
using namespace std;

void Init(int n, int* arr)//初始化
{
	for (int i = 0; i < n; i++) {
		arr[i] = i + 1;
	}
}


int simple(int n, int* arr) {//平凡算法
	int sum = 0;
	for (int i = 0; i < n; i++) {
		sum += arr[i];
	}
	return sum;
}


int duolu(int n, int* a)//优化算法—多链路
{
	int sum1 = 0;
	int sum2 = 0;
	int sum = 0;
	for (int i = 0; i < n; i += 2) {
		sum1 += a[i];
		sum2 += a[i + 1];

	}
	sum = sum1 + sum2;
	return sum;
}


void recursion(int n, int* b)//优化算法—递归
{
	if (n == 1)return;
	else {
		for (int i = 0; i < (n / 2); i++)
		{
			b[i] += b[n - i - 1];
		}
		n = (n + 1) / 2;
		recursion(n, b);
	}
}


void cycle(int n, int* a)//优化算法—二重循环
{
	int m;
	for (int m = n; m > 1; ++m /= 2) //++m避免m为奇数时导致最末端一个数丢失
	{
		for (int i = 0; i < m / 2; i++)
			a[i] = a[i + 2] + a[i * 2 + 1];
		if (m % 4 == 1)a[m / 4] = a[m - 1];
	}
}


int unroll_duolu(int n, int* a)
{
	int sum1 = 0, sum2 = 0;
	for (int i = 0; i < n; i += 4) {
		sum1 += a[i] + a[i + 1];
		sum2 += a[i + 2] + a[i + 3];

	}
	return sum1 + sum2;
}

void unroll_cycle(int n, int* a)
{
	for (int m = n; m > 1; ++m /= 2) {
		for (int i = 0; i < m / 2; i++) {
			a[i] = a[i * 2] + a[i * 2 + 1];
			a[i + 1] = a[i * 2 + 2] + a[i * 2 + 3];

		}
	}
}


void unroll_recursion(int n, int* b)
{
	if (n == 1)return;
	else {
		for (int i = 0; i < (n / 2); i += 2)
		{
			b[i] += b[n - i - 1];
			b[i + 1] += b[n - i - 2];

		}
		n = (n + 1) / 2;
		unroll_recursion(n, b);
	}
}

int main()
{


		for (int n = 100; n <= 204800*32; n = n * 2) {

			int* arr = new int[n];
			Init(n, arr);

			LARGE_INTEGER litmp;
			LONGLONG QPart1, QPart2;
			double dfMinus, dfFreq, dfTim;
			QueryPerformanceFrequency(&litmp);
			dfFreq = (double)litmp.QuadPart;// 获得计数器的时钟频率
			QueryPerformanceCounter(&litmp);
			QPart1 = litmp.QuadPart;// 获得初始值

			for (int i = 0; i < 40; i++) {
				//int sum=simple(n, arr);
                //int sum=duolu(n, arr);
                //recursion(n, arr);
                //cycle(n, arr);
				//int sum = unroll_duolu(n, arr);
				//unroll_cycle(n, arr);
				unroll_recursion(n,arr);
				QueryPerformanceCounter(&litmp);
				QPart2 = litmp.QuadPart;//获得中止值
				dfMinus = (double)(QPart2 - QPart1);
				dfTim = dfMinus * 1000 / dfFreq;
				//cout << dfTim << endl;// 获得对应的时间值，单位为秒
			}
			cout << dfTim / 40 << endl;
		}

	return 0;
}