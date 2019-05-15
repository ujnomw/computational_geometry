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
	Node& pop()
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
    Node * getBase()
    {
        return base;
    }
    friend class Polygon;

};
  #endif // TEMPL_STACK_H
