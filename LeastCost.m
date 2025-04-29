% Least cost method
% taking inputs
supply = input("Enter the supply vector: "); % enter as column
demand = input("Enter the demand vector: "); % enter as row
cost = input("Enter the cost matrix: ");
%cost_temp = cost;
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
        % supply is more than demand
        % then we add one more demand column 
        temp = zeros(1,1);
        temp(1,1) = sum1-sum2;
        demand = [demand,temp];
        temp2 = zeros(m,1);
        cost = [cost,temp2];
    else
        % demand is more than supply
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
cost_temp = cost;
% now we do the least cost method part as this is balanced now
% if both row and coln are finished we can cross any of the row
allocations = zeros(size(cost)); 
supply1 = supply;
demand1 = demand;

while any(supply1>0) && any(demand1>0)
    [min_val,idx] = min(cost(:)); % finding min cost in the matrix
    [row,col] = ind2sub(size(cost),idx); % convert index to row,col
    alloc_amt = min(supply1(row),demand1(col));
    allocations(row,col) = alloc_amt;
    supply1(row) = supply1(row) - alloc_amt;
    demand1(col) = demand1(col)-alloc_amt;

    % 3 cases of crossing the row/col-->handeled in two cases only
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
