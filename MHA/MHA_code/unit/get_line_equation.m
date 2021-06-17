function L =  get_line_equation(X1,Y1,X2,Y2)
a = Y2 - Y1;
b = X1 - X2;
c = X2*Y1 - X1*Y2;
L = [a,b,c];
end