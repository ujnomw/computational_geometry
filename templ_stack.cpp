/*//#include "templ_stack.h"
template <typename Node >
void Stack<Node>::push(Node value)
{
	base[size] = value;
	size++;
}
template <typename Node >
Stack<Node>::Stack(Node * a)
{
	base = a;
	size = 0;
}
template <typename Node >
Node Stack<Node>::pop()
{
	size--;
	return base[size];
}
template <typename Node >
Node Stack<Node>::gettop()
{
	return base[size - 1];
}
template <typename Node >
Node Stack<Node>::NextToTop()
{
  return base[size -2];
}*/

