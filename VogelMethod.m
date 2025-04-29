% Vogel Approximation method
% taking inputs
supply = input("Enter the supply vector: "); % enter as column
demand = input("Enter the demand vector: "); % enter as row
cost = input("Enter the cost matrix: ");
% first checking if this is balanced or imbalanced
sum1 = sum(supply);
sum2 = sum(demand);
disp("supply and demand comparision: ");
disp(sum1);
disp(sum2);
[m,n] = size(cost);
if sum(supply)~=sum(demand)
    disp("This is not balanced!");
    % we need to make this balanced and as no priority given so we add
    % dummy column with all 0 cost
    if sum1>sum2
        % supply more than demand
        % then we add one more demand column 
        temp = zeros(1,1);
        temp(1,1) = sum1-sum2;
        demand = [demand,temp];
        temp2 = zeros(m,1);
        cost = [cost,temp2];
    else
        % demand more than supply
        % then add one more supply column
        temp = zeros(1,1);
        temp(1,1) = sum2-sum1;
        supply = [supply;temp];
        temp2 = zeros(1,n);
        cost = [cost;temp2];
    end
end
disp("Updated Supply:");
disp(supply);
disp("Updated Demand:");
disp(demand);
disp("Updated Cost Matrix:");
disp(cost);
cost_temp = cost; % for using later on for multiplication
% new dimensions
[m,n] = size(cost);
allocations = zeros(size(cost)); 
supply1 = supply;
demand1 = demand;

while any(supply1>0) && any(demand1>0)
    % first we caculate the row and coln penalties 
    % and then choose the max and get the min cost cell corresponding to
    % that row or col
    % and there we do the allocaion
    % and put the cost to inf
    % and wherever tie in the penalty are there choose row only
    % whenever tie is there in the min cost cell choose any one 

    row_penalty = zeros(m,1);
    col_penalty = zeros(1,n);
    for i = 1 : m
        valid_cost = cost(i,demand1>0); % only consider the active columns
        % we could also have simply ignore the inf as after allocating we
        % put inf there
        if length(valid_cost)>1
            sorted_costs = sort(valid_cost);
            row_penalty(i) = sorted_costs(2) - sorted_costs(1);
        else
            row_penalty(i) = 0;
        end
    end
    for j = 1:n
        valid_cost = cost(supply1>0,j);
        if length(valid_cost)>1
            sorted_costs = sort(valid_cost);
            col_penalty(j) = sorted_costs(2) - sorted_costs(1);
        else
            col_penalty(j) = 0;
        end
    end
    % finding the max penalty
    [max_row_penalty,row_idx] = max(row_penalty);
    [max_col_penalty,col_idx] = max(col_penalty);
    if max_row_penalty>=max_col_penalty
        % find the min cost in the selected row
        [~,min_col] = min(cost(row_idx,:));
        i = row_idx;
        j = min_col;
    else
        [~,min_row] = min(cost(:,col_idx));
        i = min_row;
        j = col_idx;
    end
        % now we will do the allocations
        alloc_amt = min(supply1(i),demand1(j));
        allocations(i,j) = alloc_amt;
        supply1(i) = supply1(i) - alloc_amt;
        demand1(j) = demand(j) - alloc_amt;
        % again handling the 3 cases using two cases (if-else)
        if supply1(row)==0
        % i choose to cut the row
        cost(row,:) = inf;
        else 
        cost(:,col) = inf; 
        end
end
% Display final allocation
disp("Final Allocation:");
disp(allocations);

% Compute total transportation cost
total_cost = sum(sum(allocations .* cost_temp));
disp("Total Transportation Cost:");
disp(total_cost);