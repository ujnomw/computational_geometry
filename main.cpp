#include <iostream>
#include "geometry.h"
#include <iomanip>
int main() {

    int n;
    std::cin >> n;
    Point* s = new Point[n];
    double x,y;
    for (int i = 0; i< n; i++)
    {
       // std::cin >> x >> y;
        scanf("%lf""%lf",&x,&y);
        Point temp(x,y);
        s[i] = temp;
    }
    Polygon poly;
    /*std::cout.setf(std::ios::fixed);
    std::setprecision(0);*/
    poly.andrewHull(s,n);
    //poly.JarvisHull(s,n);
    ull sq = poly.Square();
    sq = llabs(sq);
    //sq /= 2;
    std::cout << (sq / 2 )<< '.' << (sq % 2) * 5 << '\n';

}
