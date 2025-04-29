% Identify artificial variable indices
artificial_indices = find(c == -1);
disp("These are the artificial variables: ")
disp(artificial_indices)

% Remove artificial variable columns from A
table(:, artificial_indices) = [];
disp("This is table without artificial variables");
display(table);

% Deploy the new objective function for Phase II
% Phase II objective function (keeping only original variables)
c2 = [C, extra_mat_for_objective_phase2];

% **New Cb1: Assign correct cost coefficients for current basic variables**
Cb1 = c2(basic_vars);
disp("Updated Cb1 for Phase II:");
disp(Cb1);

% **Extract A matrix from table (excluding RHS column)**
A = table(:, 1:end-1); % Excludes the last column (b)

% **Compute Initial Zj - Cj for Phase II**
Zj_minus_Cj = Cb1 * A - c2;
disp("Initial Zj - Cj for Phase II:");
disp(Zj_minus_Cj);

% **Simplex Iterations for Phase II**
no_of_iter = 0;
max_iterations = nchoosek(no_of_var + no_of_constraints, no_of_constraints);

while any(Zj_minus_Cj < 0)  % Continue until all Zj - Cj ≥ 0
    no_of_iter = no_of_iter + 1;
    fprintf("\nIteration %d:\n", no_of_iter);

    % Prevent infinite loops
    if no_of_iter > max_iterations
        disp("No solution exists.");
        return;
    end

    % Find entering variable (most negative Zj - Cj)
    [most_negative, entering_col] = min(Zj_minus_Cj);

    % Compute pivot row using the minimum positive ratio test
    ratios = inf(no_of_constraints, 1); % Initialize with infinity
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

    % Find the pivot row (minimum ratio)
    [min_ratio, pivot_row] = min(ratios);

    % Update the basic variable
    fprintf("Pivot element is at row %d, column %d.\n", pivot_row, entering_col);
    basic_vars(pivot_row) = entering_col;  % Update basic variables set

    % Normalize the pivot row
    pivot_element = table(pivot_row, entering_col);
    table(pivot_row, :) = table(pivot_row, :) / pivot_element;

    % Update other rows to make the pivot column zero
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table(i, :) = table(i, :) - table(i, entering_col) * table(pivot_row, :);
        end
    end

    % **Update Cb1 with new basic variable coefficients**
    Cb1(pivot_row) = c2(entering_col);

    % **Recalculate Zj - Cj using the updated basic variables**
    Zj_minus_Cj = Cb1 * table(:, 1:end-1) - c2;

    % Display updated table and Zj - Cj
    disp("Updated Simplex Table:");
    disp(table);
    disp("Zj - Cj:");
    disp(Zj_minus_Cj);
end

% **Extract the Optimal Solution**
disp("Optimal Solution:");
X = zeros(1, no_of_var + no_of_constraints); % Initialize solution vector
for i = 1:no_of_constraints
    basic_var_index = basic_vars(i);
    X(basic_var_index) = table(i, end); % Assign solution value
end

disp("Decision Variables:");
disp(X(1:no_of_var)); % Display only original variables
disp("Optimal Objective Value:");
disp(Cb1 * table(:, end)); % Compute optimal objective function value
