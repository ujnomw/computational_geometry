#ifndef TEMPL_STACK_H
#define TEMPL_STACK_H
template <typename Node>
class Stack
{
private:
	int size;
	Node * base;
public:
	void push(Node value)
	{
		base[size] = value;
		size++;
	}
	Stack(Node * a)
	{
		base = a;
		size = 0;
	}
	Node pop()
	{
		size--;
		return base[size];
	}
	Node gettop()
	{
		return base[size - 1];
	}
	Node NextToTop()
	{
	  return base[size -2];
	}
	int getsize()
	{
	  return size;
	}

};
  #endif // TEMPL_STACK_H
typedef unsigned long long ull;
  class Vector;
  class Point;
  class Segment;
  class Line;
  class Ray;
  class Vector
  {
  private:
    double x;
    double y;
  public:
      double getx();
      double gety();
      Vector();
      Vector(double f,double d);
      //Vector( Vector&);
      Vector(Point a, Point b);
      Vector(Point g);
      Vector operator+( const Vector&  v1);
      double operator*( Vector&  v1); //скалярное произвадение
      Vector operator%(double value); // умножение на число
      Vector operator-();
      double length();
      double PseudoscalarProduct(Vector&  v1); //псевдоскалярное произведние
      friend class Point;
      friend class Segment;
      friend class Line;
  };
  class Shape
  {
  public:
    virtual void shift(Vector&  v1) = 0;
    virtual bool itContains(Point p1) = 0;
    virtual bool CrossingSegment(Segment s1) = 0;
  };
  class Point: public Shape
  {
  private:
    double x;
    double y;
  public:
    Point() {x = 0; y =0;}
    double getx();
    double gety();
    Point& operator=(Point&);
    bool operator<(Point);
    Point(double, double);
    void shift(Vector&  v1);
    bool itContains(Point p1);
    bool CrossingSegment(Segment s1);
    double distanceToLine(Line l1);
    double distanceToRay( Ray r);
    double distanceToSegment(Segment s);
    friend class Segment;
    friend class Vector;
    friend class Polygon;

  };
  class Segment: public Shape
  {
  private:
    Point p1;
    Point p2;
  public:
    Point getp1() { return p1;}
    Point getp2() {return p2;}
    Segment(Point,Point);
    void shift(Vector&  v1);
    bool itContains(Point p1);
    bool CrossingSegment(Segment s1);
    double distanceToSegmnet(Segment s);
  };
  class Line: public Shape
  {
  private:
    double a;
    double b;
    double c;
    Vector direction;
  public:
    Point P_of_Cross(Line l1);
    Vector getDir() {return direction;}
    double getDirX(){ return direction.getx();}
    double getDirY(){return direction.gety();}
    Line(double , double , double );
    Line(Point, Point);
    void shift(Vector&  v1);
    bool itContains(Point p1);
    bool CrossingSegment(Segment s1);
    bool CrossingLine(Line l1);
    double getA(){return a;}
    double getB(){return b;}
    double getC(){return c;}

  };
  class Ray: public Shape
  {
  private:
      Vector dir;
      Point p;
  public:
      Ray(Point, Vector);
      void shift(Vector&  v1);
      bool itContains(Point p1);
      bool CrossingSegment(Segment s1);
      Vector getdir(){return dir;}
      Point getPoint(){return p;}

  };
  class Polygon: public Shape
  {
    unsigned long long numberOfPoints;
    Point* base;
  public:
    Polygon(unsigned long long);
    Polygon() {};
      void shift(Vector&  v1);
    bool itContains(Point p1);
    bool CrossingSegment(Segment s1);
    bool Clockwise();
    Point getPoint(ull i);
    ull getNumber() {return numberOfPoints;}
    bool isConvex();
    double Square();
    Point* convexHall();
    void JarvisHull(Point * s, int n);
  };















//#include "templ_stack.h"
//#include "geometry.h"
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;
//myfuncs
double dabs(double a)
{
  if (a < 0)
    a=a*(-1.);
  return a;
}

//VECTOR
double Vector::getx()
{
  return x;
}
double Vector::gety()
{
  return y;
}
Vector::Vector(Point a)
{
  x=a.getx();
  y=a.gety();
}
Vector::Vector()
{
    x = 0;
    y = 0;
}
Vector::Vector(double a,double b)
{
    x=a;
    y=b;
}

//Vector( Vector&)
Vector Vector::operator-()
{
  Vector v(-x,-y);
  return v;
}

Vector Vector::operator+( const Vector& v1)
{
    Vector vf((v1.x+x), (v1.y + y));
    return vf;
}
double Vector::operator*( Vector&  v1) //скалярное произвадение
{
  return v1.getx()*x + v1.gety()*y;
}
Vector Vector::operator%( double value) // умножение на число
{
  Vector vf(x*value, y*value);
  return vf;
}
double Vector::PseudoscalarProduct(Vector& v1) //псевдоскалярное произведние
{
  return (x*v1.gety()-y*v1.getx());
}
double Vector::length()
{
  return std::sqrt(x*x +y*y);
}
Vector::Vector(Point head, Point tail)
{
  x=tail.getx() - head.getx();
  y=tail.gety() - head.gety();
}
/*#######################################################################################################*/

//POINT
Point& Point::operator=(Point& a)
{
  x=a.getx();
  y=a.gety();
  return *this;
}

Point::Point(double a, double b)
{
    x=a;
    y=b;
}
double Point::getx()
{
    return x;
}
double Point::gety()
{
    return y;
}
void Point::shift(Vector&  v1)
{
  x=x+v1.getx();
  y=y+v1.gety();
}
bool Point::itContains(Point p1)
{
  return ((p1.getx() == x) and (p1.gety() == y));
}
bool Point::CrossingSegment(Segment s1)
{
    Vector v1(s1.getp1(), s1.getp2());
    Vector v2(s1.getp1(),*this);
    Vector v3(*this, s1.getp1());
    Vector v4(*this, s1.getp2());
    return ((v1.PseudoscalarProduct(v2) == 0) && (v3*v4 <= 0));
}
double Point::distanceToLine(Line l1)
{
    return std::abs((l1.getA()*x+l1.getB()*y +l1.getC()))/ sqrt(l1.getA()*l1.getA() + l1.getB()*l1.getB());
}
double Point::distanceToRay(Ray r)
{
    Vector v1(r.getPoint(), *this );
    Vector v2 = r.getdir();
    if (v1*v2 >= 0)
    {
        Point p1 = r.getPoint();
        p1.shift(v2);
        Line l1(p1, r.getPoint());
        return this->distanceToLine(l1);
    }
    else
        return v1.length();
}
double Point::distanceToSegment(Segment s)
{
  Vector v1(s.getp1(),*this);
  Vector v2(s.getp1(), s.getp2());
  Vector v3(s.getp2(), *this);
  Vector v4 = -v2;
  if ((v1*v2 < 0) || (v3*v4 < 0))
    {
      return min<double> (v1.length(), v3.length());
    }
  else
    {
      Line l1(s.getp1(),s.getp2());
      return this->distanceToLine(l1);
    }
}
bool Point::operator<( Point b)
{
  return ((x<b.getx()) && (y<b.gety()));
}

/*#######################################################################################################*/
//SEGMENT
Segment::Segment(Point a,Point b): p1(a), p2(b) {}
void Segment::shift(Vector& v1)
{
  p1.shift( v1);
  p2.shift( v1);
}
bool Segment::itContains(Point p1)
{
  return p1.CrossingSegment(*this);
}
bool Segment::CrossingSegment(Segment s1)
{
  Vector v1(s1.getp1(),s1.getp2()), v2(s1.getp1(),p2), v3(s1.getp1(),s1.getp2()), v4(s1.getp1(),p1), v5(p1,p2), v6(p1,s1.getp1()), v7(p1,p2), v8(p1,s1.getp2());
  double psp1 = v1.PseudoscalarProduct(v2);
  double psp2 = v3.PseudoscalarProduct(v4);
  double psp3 = v5.PseudoscalarProduct(v6);
  double psp4 = v7.PseudoscalarProduct(v8);
  return (((psp1*psp2 < 0) && (psp3*psp4 < 0)) || (s1.itContains(p1)) || (s1.itContains(p2)) || (this->itContains(s1.getp1())) || (this->itContains(s1.getp2())));
}
double Segment::distanceToSegmnet(Segment s)
{
  if (s.CrossingSegment(*this))
    return 0;
  else
    {
      return min(min(p1.distanceToSegment(s),p2.distanceToSegment(s)),min(s.getp1().distanceToSegment(*this),s.getp2().distanceToSegment(*this)));
    }
}

/*#######################################################################################################*/

//LINE
Line::Line(double a1, double b1, double c1)
{
    a = a1;
    b = b1;
    c = c1;
   Vector v1((-1)*b1, a1);
    direction = v1;
}
Line::Line(Point a, Point b)
{
   this->a = a.gety()-b.gety();
   this->b = b.getx()-a.getx();
   this->c = a.getx()*b.gety() - a.gety()*b.getx();
}

void Line::shift(Vector& v1)
{
  c = c - a*v1.getx() - b*v1.gety();
}
bool Line::itContains(Point p1)
{
  return (a*p1.getx() + b*p1.gety() + c == 0);
}
bool Line::CrossingSegment(Segment s1)
{
    Vector v1(s1.getp1());
    Vector v2(s1.getp2());
    return ((direction.PseudoscalarProduct(v1)*direction.PseudoscalarProduct(v2))<0);
}
bool Line::CrossingLine(Line l1)
{
    Vector v1=l1.getDir();
    return (this->getDir().PseudoscalarProduct(v1) != 0);
}
Point Line::P_of_Cross(Line l1)
{
    double a1 = a;
    double b1 = b;
    double c1 = c;
    double a2 = l1.getA();
    double b2 = l1.getB();
    double c2 = l1.getC();
    double x = -(c1*b2 - c2*b1)/(a1*b2 - a2*b1);
    double y = -(a1*c2 - a2*c1)/(a1*b2 - a2*b1);
    Point p(x,y);
    return p;
}

/*#######################################################################################################*/

//RAY
Ray::Ray(Point p1, Vector v1)
{
    this->dir = v1;
    this->p = p1;
}
void Ray::shift(Vector&  v1)
{
    p.shift(v1);
}
bool Ray::itContains(Point p1)
{
    Vector v1(p,p1);
    return ((p1.itContains(p)) || (((dir.PseudoscalarProduct(v1))==0) && (dir*v1 >= 0)) );
}

bool Ray::CrossingSegment(Segment s1)
{
    Vector v1(p,s1.getp1());
    Vector v2(p,s1.getp2());
    Vector vdown,vup;
    if (v1.PseudoscalarProduct(dir) > 0)
    {
        vdown = v1;
        vup = v2;
    }
    else
    {
        vdown = v2;
        vup = v1;
    }

    return ( (((dir.PseudoscalarProduct(v1)*dir.PseudoscalarProduct(v2))<0) && (vdown.PseudoscalarProduct(vup)>0)) || (this->itContains(s1.getp1())) || (this->itContains(s1.getp2())) || (s1.itContains(p)) );
}

/*#######################################################################################################*/
//POLYGONE
Polygon::Polygon(ull a)
{
  int x1, y1;
  base = new Point[a];
  numberOfPoints = a;
  for (ull i=0; i< a; i++)
    {
      //std::cout << "here";
      std::cin >> x1 >> y1;
      Point cur(x1,y1);
      base[i]=cur;
    }
}
void Polygon::shift(Vector&  v1)
{
  for (ull i = 0; i < numberOfPoints; i++ )
    {
      base[i].shift(v1);
    }
}

bool Polygon::itContains(Point p1)
{
  Vector v(0,1);
  Ray r(p1, v);
  for (ull i = 0; i < numberOfPoints; i++ )
    {
      if (r.itContains(base[i]))
        {
          double d = 1/1001.;
          base[i].x -= d;
          base[i].y -= d;
        }
    }
  int crossingNumber = 0;
 for (ull i = 0; i < numberOfPoints;i++)
   {
     Segment s(base[i],base[i % numberOfPoints]);
     if (r.CrossingSegment(s))
       crossingNumber++;
   }
 if (crossingNumber % 2 == 0)
   return true;
 return false;
}

bool Polygon::CrossingSegment(Segment s1)
{
  for (ull i = 0; i < numberOfPoints-1; i++ )
    {
      Segment s(base[i],base[i+1]);
      if (s1.CrossingSegment(s))
        return true;
    }
  return false;
}

bool Polygon::Clockwise()
{
  Vector v1(this->getPoint(0), this->getPoint(1));
  Vector v2(this->getPoint(0), this->getPoint(this->getNumber()-1));

  if (v1.PseudoscalarProduct(v2) < 0)
    return true;
  else
    return false;
}
Point Polygon::getPoint(ull i)
{
  return base[i];
}
bool Polygon::isConvex()
{
  double mult;
  ull n = getNumber() -2;
  int *product=new int[numberOfPoints];
  for (ull i = 0; i<n; i++)
    {
      Vector v1(this->getPoint(i), this->getPoint(i+1));
      Vector v2(this->getPoint(i), this->getPoint(i+2));
      mult = v1.PseudoscalarProduct(v2);
      if (mult != 0)
        product[i] = int(mult/dabs(mult));
      else
        product[i]=0;
    }
  Vector v1(this->getPoint(numberOfPoints-2), this->getPoint(numberOfPoints-1));
  Vector v2(this->getPoint(numberOfPoints-1), this->getPoint(0));
  Vector v3(this->getPoint(0),this->getPoint(1));
  mult = v1.PseudoscalarProduct(v2);
  if (mult != 0)
    product[n] = int(mult/dabs(mult));
  else
    product[n]=0;
  mult = v2.PseudoscalarProduct(v3);
  if (mult != 0)
    product[n+1] = int(mult/dabs(mult));
  else
    product[n+1]=0;
  int samepo = 0, samepr =0;
  for (ull i = 0; i<numberOfPoints;i++)
    if (product[i] == 0)
      {
        samepo++;
        samepr++;
      }
  else
      {
        if (product[i] > 0)
          samepr++;
        else
          samepo++;
      }
 if ((samepo == numberOfPoints) || (samepr == numberOfPoints))
   return true;
 return false;
}
double Polygon::Square()
{
  double sum;
  for (ull i = 0;i < numberOfPoints;i++)
    {
      Vector v1(this->getPoint(i));
      Vector v2(this->getPoint((i+1) % numberOfPoints));
      sum += v1.PseudoscalarProduct(v2);
      //std::cout << v1.PseudoscalarProduct(v2) << '\n';
    }
  return abs(sum)/2.;
}
Point* Polygon::convexHall()
{
  ull n = numberOfPoints;
  for (ull i = 0; i < n; i++)
    {
      if (base[i] < base[0])
        {
          Point temp = base[0];
          base[0] = base[i];
          base[i] = temp;
        }
    }
  for (ull i = 1; i < n; i++)
    {
      ull j = i;
      while (j >= 1)
        {
          //std::cout << "woo hoo" << '\n';
          Vector v1(this->getPoint(0),this->getPoint(j-1));
          Vector v2(this->getPoint(j-1),this->getPoint(j));
          if (v1.PseudoscalarProduct(v2) < 0)
            {
              Point temp = base[j];
              base[j] = base[j-1];
              base[j-1] = temp;
              break;
            }
          j--;
        }
    }
  Point *Hull;
  Hull = new Point[numberOfPoints];
  Stack <Point> Route(Hull);
  for (ull i = 0; i < numberOfPoints; i++)
    {
      std::cout << this->getPoint(i).getx() << ' ' << this->getPoint(i).gety() << '\n';
    }
  Route.push(base[0]);
  Route.push(base[1]);
  for (ull i = 1; i < numberOfPoints; i++)
    {
      Vector v1(Route.NextToTop(),Route.gettop());
      Vector v2(Route.gettop(),this->getPoint(i));
      double mult = v1.PseudoscalarProduct(v2);
      while (mult < 0)
        {
          Vector v1(Route.NextToTop(),Route.gettop());
          Vector v2(Route.gettop(),this->getPoint(i));
          mult = v1.PseudoscalarProduct(v2);
          Route.pop();
        }
      Route.push(this->getPoint(i));
    }
  n = Route.getsize();
  for (ull i = 0;i < n ; i++)
    {
     Point p = Route.pop();
     std::cout << p.getx() << ' ' << p.gety() << '\n';
    }
}


void Polygon::JarvisHull(Point *s, int n)
{
    int pp = 0;
    Point p0= s[pp]; // lowest and the most right point
    for (int i = 1; i < n; i++)
    {
        if (s[i].y < p0.y)
        {
            pp = i;
            p0 = s[pp];
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (s[i].y == p0.y && s[i].x > p0.x)
        {
            pp = i;
            p0 = s[pp];
        }
    }
    Point pi = p0;
    /*Point temp = s[0];
    s[0] = s[pp];
    s[pp] = temp;*/
    int k = -1;
    //Point * result = new Point[n];
    //result[0] = p0;
    do
    {
        k++;
        for (int i = k; i < n; i++)
        {
            Vector v2(pi,s[i]);
            Vector v1(pi,s[k]);
            if (v2.PseudoscalarProduct(v1) > 0 || (v2.PseudoscalarProduct(v1) == 0 && v2.length() > v1.length()))
            {
                Point temp = s[i];
                s[i] = s[k];
                s[k] = temp;
            }
        }
        pi = s[k];
    }
    while  ( pi.x != p0.x || pi.y != p0.y );
    Point * res = new Point[k+1];
    //res[0]= p0;
    for (int i =0 ; i < k+1; i++)
        res[i] = s[i];
    base = res;
    std::cout << k+1 << '\n';
    for (int i = 0; i < k+1; i++)
        std::cout << s[i].x << ' ' << s[i].y << '\n';
    //delete [] s;
    numberOfPoints = k+1;

}






#include <iostream>
//#include "geometry.h"
int main() {
    int n;
    std::cin >> n;
    Point* s = new Point[n];
    int x,y;
    for (int i = 0; i< n; i++)
    {
        std::cin >> x >> y;
        Point temp(x,y);
        s[i] = temp;
    }
    Polygon poly;
    std::cout.precision(0);
    poly.JarvisHull(s,n);
   std::cout.precision(1000000);
    std::cout << poly.Square() << '\n';

}
