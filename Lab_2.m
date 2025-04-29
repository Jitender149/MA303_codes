% Lab 2
% linprog solve minimisation problem only
% so convert accordingly
% a)
% taking input
A = input("Enter the matrix A");
b = input("Enter the matrix b");
c = input("Enter the matrix C");
% c to be input as column vector 

[X,z] = linprog(c,A,b);
disp(X);
disp(z); % see how to adjust -ve sign
[row,col] = size(b);
hold on
x = 0:0.1:10;
for i = 1: constraints
    y = (b(i)-A(i,1)*x)/A(i,2);
    plot(x,y);
end


% how to plot the region --- shaded region to be given by the matlab
% what if there is lower bound condition on the decision problem