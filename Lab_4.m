% Lab 4 -> Assumption: All RHS are positive, and all constraints are "<="

% Input
no_of_var = input("Enter the number of variables: ");
no_of_constraints = input("Enter the number of constraints: ");
A = input("Enter the matrix A: ");
b = input("Enter the constant matrix (RHS): ");
c = input("Enter the coefficients of the objective function: ");

% Add Slack Variables
iden_mat = eye(no_of_constraints); % Identity matrix for slack variables
A = [A, iden_mat]; % Append slack variables to A

% Extend the objective function to include coefficients for slack variables (0s)
c = [c, zeros(1, no_of_constraints)];

% Initialize Simplex Table
table = [A, b]; % Combine A and b to form the table
Cb = zeros(1, no_of_constraints); % Coefficients of the basic variables (initially slack variables)
Zj_minus_Cj = Cb * A - c; % Initial Zj - Cj row, remember this is row !!!
disp("Initial Simplex Table:");
disp(table);

% Initialize iteration counter
no_of_iter = 0;

% Start Simplex Iterations
while any(Zj_minus_Cj < 0) % Continue until all Zj - Cj >= 0 (optimality condition)
    no_of_iter = no_of_iter + 1;
    fprintf("\nIteration %d:\n", no_of_iter);
    if no_of_iter> nchoosek(no_of_var+no_of_constraints,no_of_constraints);
        disp("No solution exist.")
        % we need to make n wrt to standard version ie after adding slack variables.
        return
    end
    % Find entering variable (most negative Zj - Cj)
    [most_negative, entering_col] = min(Zj_minus_Cj);

    % Compute ratios for the pivot row
    ratios = zeros(no_of_constraints, 1); % making separate matrix for this so that we can operate nicely
    for i = 1:no_of_constraints
        if table(i, entering_col) > 0
            ratios(i) = table(i, end) / table(i, entering_col);
        else
            ratios(i) = inf; % Ignore negative or zero pivot column entries
        end
    end

    % Check for unbounded solution
    if all(ratios == inf)
        disp("The solution is unbounded.");
        return;
    end

    % Find the pivot row (minimum ratio)
    [min_ratio, pivot_row] = min(ratios);

    % Update the basic variable
    fprintf("Pivot element is at row %d, column %d.\n", pivot_row, entering_col);

    % Normalize the pivot row
    pivot_element = table(pivot_row, entering_col);
    table(pivot_row, :) = table(pivot_row, :) / pivot_element; %refers to the entire pivot row (all columns of the row specified by pivot_row

    % Update other rows to make the pivot column zero
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table(i, :) = table(i, :) - table(i, entering_col) * table(pivot_row, :);
        end
    end

    % Update the coefficients of basic variables (Cb)
    Cb(pivot_row) = c(entering_col);

    % Recalculate Zj - Cj
    Zj_minus_Cj = Cb * table(:, 1:end-1) - c;

    % Display the updated table and Zj - Cj
    disp("Updated Simplex Table:");
    disp(table);
    disp("Zj - Cj:");
    disp(Zj_minus_Cj);
end

% Display the final solution ---> SEE THIS ONCE
disp("Optimal Solution:");
X = zeros(1, no_of_var + no_of_constraints); % Initialize solution vector
for i = 1:no_of_constraints
    basic_var_index = find(table(i, 1:no_of_var + no_of_constraints) == 1);
    % ofc we want to skip the last column ie RHS --->IMP 
    if ~isempty(basic_var_index)
        X(basic_var_index) = table(i, end); % Assign solution values
    end
end

disp("Decision Variables:");
disp(X(1:no_of_var)); % Display only original variables
disp("Optimal Objective Value:");
disp(Cb * table(:, end)); % Objective function value
