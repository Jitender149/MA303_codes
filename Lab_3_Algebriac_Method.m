% Lab 2 -> Algebraic Method
% Make a generalized code for solving feasible solutions using basic matrices

% Inputs
A = input("Enter the matrix A: "); % Constraint coefficient matrix
b = input("Enter the matrix B: "); % Right-hand side vector
c = input("Enter the coefficients of the objective function: "); % Objective function coefficients

% Get dimensions of A
[row, col] = size(A); 
m = row; % Number of constraints
n = col; % Number of variables

% Generate all combinations of m columns from n
combinations = nchoosek(1:n, m); 
disp("All possible combinations of columns:");
disp(combinations);

% Initialize matrix to store feasible solutions
feasSolMat = []; 

% Iterate through all combinations to generate submatrices and basic solutions
for i = 1:size(combinations, 1)
    selectCol = combinations(i, :);  % Current combination of columns    WHAT IS THE MEANING OF THIS :
    subMat = A(:, selectCol);        % Submatrix for this combination
    disp("This is the basic matrix:");
    disp(subMat);
    
    % Check if subMat is non-singular
    if det(subMat) ~= 0
        basicSol = subMat \ b; % Calculate basic solution
        
        % Check feasibility: all elements of the solution should be >= 0
        if all(basicSol >= 0)
            % Create the solution vector wrt all variables
            x = zeros(n, 1);
            for temp = 1:m
                x(selectCol(temp)) = basicSol(temp);
            end
            
            % Append the feasible solution
            feasSolMat = [feasSolMat; x'];  % dont forget to make it column vector
        end
    end
end

% Display the feasible solutions
disp("Feasible solutions:");
disp(feasSolMat);

% Now we will calculte the max value of the objective function
maxi_mat = feasSolMat*(c);
ans = max(maxi_mat)

% NOTES :
% size(combinations,1) give no of rows
% and rows(combinations) give the vector of [row,col]
% combinations(i, :) selects the i-th row from the combinations matrix.
% : means "select all rows."
% selectCol specifies which columns to take.


% % Logical vector indicating selected columns
%columnMask = [true, false, true];

% Select columns 1 and 3
%E = A(:, columnMask);
