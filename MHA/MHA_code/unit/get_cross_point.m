function [X,Y] = get_cross_point(L1,L2)

a1=L1(1);b1=L1(2);c1=L1(3);
a2=L2(1);b2=L2(2);c2=L2(3);

syms x y
L1 = [num2str(a1),'*x+',num2str(b1),'*y+',num2str(c1),'=0'];
L2 = [num2str(a2),'*x+',num2str(b2),'*y+',num2str(c2),'=0'];

[X,Y] = solve(L1,L2);
X = double(X);
Y = double(Y);

end