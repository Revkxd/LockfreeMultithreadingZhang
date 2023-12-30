#include<stdio.h>

int max_of_two(int a, int b) {
    return a > b ? a : b;
}

int max_of_three(int a, int b, int c) {
    int temp = a > b ? a : b;
    return temp > c ? temp : c;
}

int max_of_four(int a, int b, int c, int d) {
    int temp1 = a > b ? a : b;
    int temp2 = c > d ? c : d;
    return temp1 > temp2 ? temp1 : temp2;
}

int max_of_five(int a, int b, int c, int d, int e) {
    int temp1 = max_of_three(a, b, c);
    int temp2 = max_of_two(d, e);
    return temp1 > temp2 ? temp1 : temp2;
}