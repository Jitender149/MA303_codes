% Lab 5 --> Big M method
clc; clear; close all;

% Input basic details
no_of_var = input("Enter the number of variables: ");
no_of_constraints = input("Enter the number of constraints: ");
LessThanEqualTo = input("Enter the no of less than equal to constraints: ");
EqualTo = input("Enter the no of equal to constraints: ");
GreaterThanEqualTo = input("Enter the no of greater than equal to constraints: ");

% Input matrices
A = input("Enter the matrix A: ");
b = input("Enter the constant matrix (RHS): ");
c = input("Enter the coefficients of the objective function: ");

% Define a large M value
M = 1e6;

% Initialize extra columns for slack, surplus, and artificial variables
extra_mat = [];
extra_mat_for_objective = zeros(1, LessThanEqualTo + EqualTo + 2 * GreaterThanEqualTo);

% Add slack variables (for <= constraints)
for i = 1 : LessThanEqualTo
    col = zeros(no_of_constraints, 1);
    col(i) = 1;
    extra_mat = [extra_mat, col];
end

% Add artificial variables (for = constraints)
for i = 1 : EqualTo
    col = zeros(no_of_constraints, 1);
    col(LessThanEqualTo + i) = 1;
    extra_mat = [extra_mat, col];
    extra_mat_for_objective(1, LessThanEqualTo + i) = -M; % Assigning -M
end

% Add surplus and artificial variables (for >= constraints)
for i = 1 : GreaterThanEqualTo
    col1 = zeros(no_of_constraints, 1);
    col2 = zeros(no_of_constraints, 1);
    col1(LessThanEqualTo + EqualTo + i) = -1; % Surplus variable
    col2(LessThanEqualTo + EqualTo + i) = 1;  % Artificial variable
    extra_mat = [extra_mat, col1, col2];
    extra_mat_for_objective(1, LessThanEqualTo + EqualTo + i + 1) = -M; % Artificial variable weight
end

% Update A and c with the extra variables
A = [A, extra_mat];
c = [c, extra_mat_for_objective];
display(c)

% Initialize Simplex Table
table = [A, b]; % Combine A and b to form the table
Cb = zeros(1, no_of_constraints);  % Initialize Cb

% Assign 0 to slack variables (for <= constraints)
Cb(1:LessThanEqualTo) = 0;

% Assign -M (e.g., -1000) to artificial variables (for = and >= constraints)
Cb(LessThanEqualTo+1:end) = -1000;  % Adjust -M value as needed

Zj_minus_Cj = Cb * A - c; % Initial Zj - Cj row

disp("Initial Simplex Table:");
disp(table);

% Simplex Iteration
no_of_iter = 0;
max_iterations = nchoosek(no_of_var + no_of_constraints, no_of_constraints);

while any(Zj_minus_Cj < 0) % Continue until all Zj - Cj >= 0
    no_of_iter = no_of_iter + 1;
    if no_of_iter > max_iterations
        disp("No solution exists.");
        return;
    end
    
    % Find entering variable (most negative Zj - Cj)
    [most_negative, entering_col] = min(Zj_minus_Cj);

    % Compute ratios for pivot row selection
    ratios = inf(no_of_constraints, 1);
    for i = 1:no_of_constraints
        if table(i, entering_col) > 0
            ratios(i) = table(i, end) / table(i, entering_col);
        end
    end
    
    % Check for unbounded solution
    if all(ratios == inf)
        disp("The solution is unbounded.");
        return;
    end

    % Find pivot row (minimum ratio)
    [min_ratio, pivot_row] = min(ratios);

    % Normalize the pivot row
    pivot_element = table(pivot_row, entering_col);
    table(pivot_row, :) = table(pivot_row, :) / pivot_element;

    % Update other rows
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table(i, :) = table(i, :) - table(i, entering_col) * table(pivot_row, :);
        end
    end

    % Update the Cb values
    Cb(pivot_row) = c(entering_col);

    % Recalculate Zj - Cj
    Zj_minus_Cj = Cb * table(:, 1:end-1) - c;

    % Display updated table
    fprintf("\nIteration %d:\n", no_of_iter);
    disp("Updated Simplex Table:");
    disp(table);
    disp("Zj - Cj:");
    disp(Zj_minus_Cj);
end

% Extract optimal solution
disp("Optimal Solution:");
X = zeros(1, no_of_var + no_of_constraints);
for i = 1:no_of_constraints
    basic_var_index = find(table(i, 1:no_of_var + no_of_constraints) == 1);
    if ~isempty(basic_var_index)
        X(basic_var_index) = table(i, end);
    end
end

% Display results
disp("Decision Variables:");
disp(X(1:no_of_var)); % Display only original variables
disp("Optimal Objective Value:");
disp(Cb * table(:, end));
