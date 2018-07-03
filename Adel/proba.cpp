//g++ proba.cpp -o proba.out -fopenmp

#include <iostream>
#include <omp.h>
#include <cmath>
#include <climits>

using namespace std;
#define PI 3.14159265358979323846264338327950288419716939937510

time_t beginning_time = 0;
time_t ending_time = 0;

void start_timing()
{
	time(&beginning_time);
}

void stop_timing()
{
	time(&ending_time);
	double difference = difftime(ending_time, beginning_time);
	tm* time_info = localtime(&beginning_time);
	cout << "Program started at:" << asctime(time_info) << endl;

	time_info = localtime(&ending_time);
	cout << "Program ended at:" << asctime(time_info) << endl;

	cout << "Program running time was " << difference << " seconds" << endl;
}

double f(double i) {
    return sqrt(1 - i * i);
}

int main() {
    start_timing();
    int i, n = INT_MAX, nsamples = 0;
    cout << n << endl;
    double dx = 1.0 / n, pi4 = 0, h = 1;
    #pragma omp parallel for reduction(+:pi4)
    for (i=0; i<n; i++) {
        double
            x = i*h,x2 = (i+1)*h,
            y = sqrt(1-x*x),y2 = sqrt(1-x2*x2),
            slope = (y-y2)/h;
        if (slope>15) slope = 15;
        int samples = 1+(int)slope, is;
        for (is=0; is<samples; is++) {
            double
                hs = h/samples,
                xs = x+ is*hs,
                ys = sqrt(1-xs*xs);
            pi4 += hs*ys;
            nsamples++;
        }
    }
    cout << "PI = " << pi4 * 4 << endl;
    if (PI/4 - pi4 < 0.000000000000000000000000000000001)
        cout << "OK" << endl;
    stop_timing();
    return 0;
}