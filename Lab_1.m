% Lab 1
A = [ 1 2 3; 4 5 6 ; 7 8 9];
B = [ 2;4;8];
disp(A);
C = [A,B]; % adding one column in the matrix
disp(C);
% adding a row 
D = [(A') B]';
disp(D);
% in matlab space and comma works same

%inverse of a matrix
E = inv(A);
disp(E);

% display only the some part of the matrix
[row,col] = size(A)
disp(A(:,1:col-1))
disp(A(1:row-1,:))
% Always check if this is invertible or not
% first check if it is square
% and the use det(A) = 0 to check

% if this is not invertible then we use Hat matrix ie (A'*A)\(A')

F = [1 2 3; 3 4 5; 5 6 7];
G = [4 5 6; 2 8 3; 10 11 12];
syms x1 x2 x3;
x = [x1; x2; x3];  % Create a symbolic column vector for x
eqn = F * x == G;   % Define the matrix equation
sol = solve(eqn, [x1, x2, x3]);  % Solve the equation
disp(sol);

% or simply do the x = ((A'*A)\(A'))*B to get answer if matrix is not
% invertible

% plotting 
x = linspace(0, 10, 100);
y = x+5;
plot(x,y);
% plotting more than one curve in 1 figure
y1 = 2*x+3;
plot(x,y);
hold on  % it holds on the previous graph
plot(x,y1);
hold off
% now make both graph side by side
% given que in class (see ss 7/1/25)
