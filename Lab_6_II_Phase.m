% II Phase Method

% Input basic details
no_of_var = input("Enter the number of variables: ");
no_of_constraints = input("Enter the number of constraints: ");
LessThanEqualTo = input("Enter the no of less than equal to constraints: ");
EqualTo = input("Enter the no of equal to constraints: ");
GreaterThanEqualTo = input("Enter the no of greater than equal to constraints: ");

% Input matrices
A = input("Enter the matrix A: ");
b = input("Enter the constant matrix (RHS): ");
C = input("Enter the coefficients of the objective function: ");
% this is captal 


% Initialize extra columns for slack, surplus, and artificial variables
extra_mat = [];
extra_mat_for_objective_phase1 = zeros(1, LessThanEqualTo + EqualTo + 2 * GreaterThanEqualTo); % we only need to assign the -1 for the artificial var
extra_mat_for_objective_phase2 = zeros(1, LessThanEqualTo + GreaterThanEqualTo);
% in phase II there will not be artificila var but ofc there can be slack,
% surplus var and their wt would be 0 only
% so its better to add only wt 0 for slack and surplus var only 
% as we would drop the artificial columns after phase I


% Add slack variables (for <= constraints)
for i = 1 : LessThanEqualTo
    col = zeros(no_of_constraints, 1);
    col(i) = 1;
    extra_mat = [extra_mat, col]; % --> this one is for the A matrix
    % objective in both phase 1 and phase 2 would have 0 only
end

% Add artificial variables (for = constraints)
for i = 1 : EqualTo
    col = zeros(no_of_constraints, 1);
    col(LessThanEqualTo + i) = 1;
    extra_mat = [extra_mat, col];
    extra_mat_for_objective_phase1(1, LessThanEqualTo + i) = -1; % Assigning -1
    % as this is = type so for phase 2 this artifical would be removed
end

% Add surplus and artificial variables (for >= constraints)
for i = 1 : GreaterThanEqualTo
    col1 = zeros(no_of_constraints, 1);
    col2 = zeros(no_of_constraints, 1);
    col1(LessThanEqualTo + EqualTo + i) = -1; % Surplus variable
    col2(LessThanEqualTo + EqualTo + i) = 1;  % Artificial variable
    extra_mat = [extra_mat, col1, col2];
    extra_mat_for_objective_phase1(1, LessThanEqualTo + EqualTo + i + 1) = -1; % Artificial variable weight
end

% Update A and c with the extra variables
A = [A, extra_mat];
c1 = zeros(1,no_of_var);
c = [c1, extra_mat_for_objective_phase1]; % as for phase I we need all zero except -1 for atrificial var
display(c)
c2 = [C, extra_mat_for_objective_phase2];
display(c2); % this one will be used in phase 2

% Track Basic Variables (Initially, these are slack & artificial variables)
basic_vars = (no_of_var + 1 : no_of_var + no_of_constraints);

% Objective function for Phase I
% only artificial var are given weights = -1 rest all put zero 
% now we need to get the basic matrix and accordingly update the A and make
% the Cb accordingly

% Initialize Simplex Table
table = [A, b]; % Combine A and b to form the table
Cb = zeros(1, size(A,1));  % Use the number of rows in A
  % Initialize Cb

% Assign 0 to slack variables (for <= constraints)
Cb(1:LessThanEqualTo) = 0;  

% Assign -1 to artificial variables (for = and ≥ constraints)
Cb(LessThanEqualTo+1 : LessThanEqualTo+EqualTo) = -1; % Artificial vars from equality constraints
Cb(LessThanEqualTo+EqualTo+1 : end) = -1; % Artificial vars from ≥ constraints

disp(Cb); % Now it should display correctly
disp(A)
Zj_minus_Cj = Cb * A - c; % Initial Zj - Cj row

disp("Initial Simplex Table:");
disp(table);

% Simplex Iteration
no_of_iter = 0;
max_iterations = nchoosek(no_of_var + no_of_constraints, no_of_constraints);

while any(Zj_minus_Cj < 0) % Continue until all Zj - Cj >= 0
    no_of_iter = no_of_iter + 1;
    if no_of_iter >= max_iterations
        disp("No solution exists.");
        return;
    end
    
    % Find entering variable (most negative Zj - Cj)
    [most_negative, entering_col] = min(Zj_minus_Cj);
       disp("Entering column is: ")
       disp(entering_col);
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
    disp("Leaving Column is: ")
    disp(basic_vars(pivot_row));  % Display the variable that will leave the basis
    % Normalize the pivot row
    pivot_element = table(pivot_row, entering_col);
    table(pivot_row, :) = table(pivot_row, :) / pivot_element;

    % Update other rows
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table(i, :) = table(i, :) - table(i, entering_col) * table(pivot_row, :);
        end
    end
    
    
    % Update Basic Variables
    basic_vars(pivot_row) = entering_col;

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
disp("These are the basic variables: ")
disp(basic_vars);
% now here phase I will end and we have to drop the artificial columns form
% A and deploy new objective values for the Phase II
% corresponding to the basic variables wrt the way in which these we
% written

% Identify artificial variable indices
artificial_indices = find(c == -1);
disp("These are the artificial variables: ")
disp(artificial_indices)
% Remove artificial variable columns from A
table(:, artificial_indices) = [];
disp("This is table without artificial variables");
display(table);
% deploy the new objective values
% C2 simply is the new objective values
% we just need to make new Cb' (coefft of basic var in phase II starting)
new_C = [C,extra_mat_for_objective_phase1]
Cb1 = new_C(basic_vars);  % isme se hmesha non artificial wale hi aeynge
disp(Cb1);
% and assuming the artificila var not present in the last table of Phase I
% so we will assign new Cb values but before that we need to know which all
% var are basic in the last table


% now again the simplex
% Extract A matrix from the table (excluding the RHS column)
A = table(:, 1:end-1); % Excludes the last column (b)

% Compute Zj - Cj for Phase II
Zj_minus_Cj = Cb1 * A - c2; % Initial Zj - Cj row

% Display the initial simplex table
disp("Initial Simplex Table:");
disp(table);


% Initialize iteration counter
no_of_iter = 0;

% Start Simplex Iterations for Phase II
while any(Zj_minus_Cj < 0) % Continue until all Zj - Cj >= 0 (optimality condition)
    no_of_iter = no_of_iter + 1;
    fprintf("\nIteration %d:\n", no_of_iter);
    
    % Prevent infinite loops (upper bound on iterations)
    if no_of_iter > nchoosek(no_of_var + no_of_constraints, no_of_constraints)
        disp("No solution exists.");
        return;
    end

    % Find entering variable (most negative Zj - Cj)
    [most_negative, entering_col] = min(Zj_minus_Cj);

    % Compute ratios for the pivot row
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
    
    % **Track Basic Variables**
    basic_vars(pivot_row) = entering_col; % Update basic variables set

    % Normalize the pivot row
    pivot_element = table(pivot_row, entering_col);
    table(pivot_row, :) = table(pivot_row, :) / pivot_element; % Normalize row

    % Update other rows to make the pivot column zero
    for i = 1:no_of_constraints
        if i ~= pivot_row
            table(i, :) = table(i, :) - table(i, entering_col) * table(pivot_row, :);
        end
    end

    % Update the coefficients of basic variables (Cb1)
    Cb1(pivot_row) = c2(entering_col);

    % Recalculate Zj - Cj using the updated basic variables
    Zj_minus_Cj = Cb1 * table(:, 1:end-1) - c2;

    % Display the updated table and Zj - Cj
    disp("Updated Simplex Table:");
    disp(table);
    disp("Zj - Cj:");
    disp(Zj_minus_Cj);
end

% **Extract the Optimal Solution**
disp("Optimal Solution:");
X = zeros(1, no_of_var + no_of_constraints); % Initialize solution vector
for i = 1:no_of_constraints
    basic_var_index = basic_vars(i); % Get basic variable index
    X(basic_var_index) = table(i, end); % Assign solution value
end

disp("Decision Variables:");
disp(X(1:no_of_var)); % Display only original variables
disp("Optimal Objective Value:");
disp(Cb1 * table(:, end)); % Compute optimal objective function value
