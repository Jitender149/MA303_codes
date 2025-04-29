% Vogel Approximation Method
% Taking inputs
supply = input("Enter the supply vector: "); % Enter as column
demand = input("Enter the demand vector: "); % Enter as row
cost = input("Enter the cost matrix: ");

% Store the original cost matrix
%cost_temp = cost;  

% Checking if the problem is balanced
sum1 = sum(supply);
sum2 = sum(demand);
disp("Supply and Demand comparison: ");
disp(sum1);
disp(sum2);

[m, n] = size(cost);

% now we will balance the matrix 
if sum1 ~= sum2
    disp("This is not balanced!");
    % we need to make this balanced and as no priority given so we add
    % dummy column with all 0 cost
    if sum1 > sum2
        % supply more than demand
        % then we add one more demand column 
        demand = [demand, sum1 - sum2];
        cost = [cost, zeros(m, 1)];
    else
        % demand more than supply
        % then add one more supply column
        supply = [supply; sum2 - sum1];
        cost = [cost; zeros(1, n)];
    end
end

disp("Updated Supply:");
disp(supply);
disp("Updated Demand:");
disp(demand);
disp("Updated Cost Matrix:");
disp(cost);

cost_temp = cost; % For final cost calculation
allocations = zeros(size(cost)); 
supply1 = supply;
demand1 = demand;

while any(supply1 > 0) && any(demand1 > 0)
    % first we caculate the row and coln penalties 
    % and then choose the max and get the min cost cell corresponding to
    % that row or col
    % and there we do the allocaion
    % and put the cost to inf
    % and wherever tie in the penalty are there choose row only
    % whenever tie is there in the min cost cell choose any o
    % Compute row and column penalties
    row_penalty = zeros(m, 1);
    col_penalty = zeros(1, n);
    for i = 1:m
        valid_costs = cost(i, demand1 > 0); % Consider active columns
        if length(valid_costs) > 1
            sorted_costs = sort(valid_costs);
            row_penalty(i) = sorted_costs(2) - sorted_costs(1);
        else
            row_penalty(i) = 0;
        end
    end
    for j = 1:n
        valid_costs = cost(supply1 > 0, j); % Consider active rows
        if length(valid_costs) > 1
            sorted_costs = sort(valid_costs);
            col_penalty(j) = sorted_costs(2) - sorted_costs(1);
        else
            col_penalty(j) = 0;
        end
    end

    % now finding the maximum penalty
    [max_row_penalty, row_idx] = max(row_penalty);
    [max_col_penalty, col_idx] = max(col_penalty);
    if max_row_penalty >= max_col_penalty
        [~, min_col] = min(cost(row_idx, :));
        i = row_idx;
        j = min_col;
    else
        [~, min_row] = min(cost(:, col_idx));
        i = min_row;
        j = col_idx;
    end

    % Now we will do the allocations of this min cell we selected
    alloc_amt = min(supply1(i), demand1(j));
    allocations(i, j) = alloc_amt;

    % Update supply and demand vector copies
    supply1(i) = supply1(i) - alloc_amt;
    demand1(j) = demand1(j) - alloc_amt;
    % again handling the 3 cases using two cases (if-else)
    if supply1(i) == 0
        cost(i, :) = inf;   
    else
        cost(:, j) = inf;
    end

    % If only one row or column is remaining ,using the least cost method
    active_rows = find(supply1 > 0);
    active_cols = find(demand1 > 0);
    if length(active_rows) == 1 || length(active_cols) == 1
        for i = active_rows
            for j = active_cols
                allocations(i, j) = min(supply1(i), demand1(j));
                supply1(i) = supply1(i) - allocations(i, j);
                demand1(j) = demand1(j) - allocations(i, j);
            end
        end
        break;
    end
end
disp("These are the final allocations: ");
disp(allocations);

% Compute final cost
total_cost = sum(sum(allocations .* cost_temp));
disp("Total Transportation Cost:");
disp(total_cost);
