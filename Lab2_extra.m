clc; clear; close all;

% Taking input
A = input("Enter the matrix A: ");
b = input("Enter the matrix b: ");
c = input("Enter the cost vector C (as a column vector): ");

% Solve the LP using linprog
[X, z] = linprog(c, A, b);

disp("Optimal Solution:");
disp(X);
disp("Optimal Objective Value:");
disp(z);

% Extract number of constraints
[row, ~] = size(A);

hold on; % Keep all plots on the same figure
x_range = 0:0.1:10; % Define x-axis range

% Plot each constraint line
for i = 1:row
    if A(i,2) ~= 0  % Avoid division by zero
        y_values = (b(i) - A(i,1) * x_range) / A(i,2); % Correct formula
        plot(x_range, y_values, 'b', 'LineWidth', 2); % Plot the constraint line
    else
        % If A(i,2) == 0, the constraint is a vertical line x = b(i)/A(i,1)
        x_fixed = b(i) / A(i,1);
        y_fixed = 0:0.1:10; % y-values for vertical line
        plot(x_fixed * ones(size(y_fixed)), y_fixed, 'r', 'LineWidth', 2);
    end
end

% Set axis limits
xlim([0 10]); 
ylim([0 10]); 

hold off;
grid on;
title("Feasible Region for Linear Programming");
xlabel("x-axis");
ylabel("y-axis");
legend("Constraint Lines");
