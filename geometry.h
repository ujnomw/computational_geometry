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
      friend bool operator<(Point& a,Point & b)
      {

      }
      long long HullCrutchPs(Vector& v1);
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
    friend int compareBeneath(const void * left, const void * right);
   /* {
        Point p1 = *(Point *)left;
        Point p2 = *(Point *)right;
        return (int)(p2.x - p1.x);
    }*/
      friend int compareAbove(const void * left, const void * right);
     /* {
          Point p1 = *(Point *)left;
          Point p2 = *(Point *)right;
          return (int)(p1.x - p2.x);
      }*/

    friend void swap(Point& a, Point& b);
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
      long long Square();
    void JarvisHull(Point * s, int n);
    void andrewHull(Point *s, int n);
      long long Perimetr();
  };















