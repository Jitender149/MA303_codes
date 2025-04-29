clc;
clear all;

% Problem Inputs
% Example: Maximize 3x + 2y
% Subject to: x + y <= 4
%             x - y >= 1   (i.e., -x + y <= -1)
%             x, y >= 0
%             x ∈ Integers

c = [3; 5];         % Objective coefficients
A = [2 4; 1 0;0 -2];    % Constraint matrix
b = [25; 8;10];        % Right-hand side
int_vars = [1,2];     % Indices of integer variables (e.g., only x)

% Convert c and b to column vectors
c = c(:);
b = b(:);

n = length(c); % Number of variables

% Root Node
root.problem.A = A;
root.problem.b = b;
root.problem.c = c;
root.problem.int_vars = int_vars;
root.problem.bounds = [zeros(n, 1), inf(n, 1)];  % Default bounds [0, ∞)
root.level = 0;

% Initialization
best_sol = [];
best_val = -inf;
queue = {root};    % Stack of nodes (for DFS)

while ~isempty(queue)
    % Pop the last node
    node = queue{end};
    queue(end) = [];

    % Solve LP relaxation at current node
    options = optimoptions('linprog', 'Display', 'none');
    [x, fval, exitflag] = linprog(-node.problem.c, node.problem.A, node.problem.b, [], [], ...
                                  node.problem.bounds(:,1), node.problem.bounds(:,2), options);
    fval = -fval;  % Convert back to maximization

    if exitflag ~= 1
        continue;  % Infeasible or unbounded
    end

    % Check if solution satisfies integrality for integer variables
    is_integer = true;
    for idx = node.problem.int_vars
        if abs(x(idx) - round(x(idx))) > 1e-5
            is_integer = false;
            break;
        end
    end

    if is_integer
        % Update best solution if better
        if fval > best_val
            best_val = fval;
            best_sol = x;
        end
        continue;
    end

    % Otherwise, branch on the most fractional variable
    frac_var = -1;
    max_frac = 0;
    for idx = node.problem.int_vars
        frac = abs(x(idx) - floor(x(idx)));
        if frac > max_frac + 1e-5 && frac < 1 - 1e-5
            max_frac = frac;
            frac_var = idx;
        end
    end

    if frac_var == -1
        continue; % No branching variable found
    end

    x_frac = x(frac_var);

    % Create two child nodes:
    % Left child: x_frac <= floor(x_frac)
    new_node1 = node;
    new_row = zeros(1, n);
    new_row(frac_var) = 1;
    new_node1.problem.A = [node.problem.A; new_row];
    new_node1.problem.b = [node.problem.b; floor(x_frac)];
    new_node1.level = node.level + 1;

    % Right child: x_frac >= ceil(x_frac)
    new_node2 = node;
    new_row = zeros(1, n);
    new_row(frac_var) = -1;
    new_node2.problem.A = [node.problem.A; new_row];
    new_node2.problem.b = [node.problem.b; -ceil(x_frac)];
    new_node2.level = node.level + 1;

    % Add children to stack
    queue{end+1} = new_node1;
    queue{end+1} = new_node2;
end

% Final Output
fprintf('\nBest Integer Solution Found:\n');
disp(best_sol);

fprintf('Optimal Integer Objective Value: %.4f\n', best_val);
