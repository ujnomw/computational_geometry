n = int(input())
lisy = []
for i in range(n):
    x = int(input())
    y = int(input())
    x /=    10000000;
    y /= 10000000;
    l1 = [x,y]
    lisy.appednd(l1)
for i in lisy:
    print(i.x, i.y, sep=' ')
