% North West corner rule
% write whats happening in the code
% handle the unbalanced case as well
a = input("Enter the supply : "); % enter as column vector
b = input("Enter the demand : "); % enter as row vector
cost = input("Enter the cost matrix : ");
n = size(a,1);
m = size(b,2);
disp("These are the no of rows in supply: ");
disp(n);
disp("These are the no of cols in demand: ");
disp(m);
% now handling the unbalanced part
tot_supply = sum(a)
tot_demand = sum(b)
if tot_demand>tot_supply
    disp("This is not balanced. ");
    % ie we need to add a supply row with all cost 0
    temp = zeros(1,1);
    temp = tot_demand-tot_supply;
    a = [a;temp];
    disp("This is new supply: ");
    disp(a);
    temp1 = zeros(1,m);
    cost = [cost;temp1];
    disp("This is new cost matrix: ");
    disp(cost);
end
if tot_supply>tot_demand
    disp("This is not balanced, ");
    % ie we need to add a demand col with all costs 0
    temp = zeros(1,1);
    temp = tot_supply-tot_demand;
    b = [b,temp];
    disp("This is new demand: ");
    disp(b);
    temp1 = zeros(1,n);
    cost = [cost,temp1];
    disp("This is new cost matrix: ");
    disp(cost);
end
[n,m] = size(cost);
t = zeros(n,m); % doing the dimesion given in the q1
i = 1;
j = 1;

while(i<=n && j<=m)
    if a(i)<b(j)
        t(i,j)=a(i);
        b(j) = b(j) - a(i);
        i = i+1;

    else
        t(i,j)= b(j);
        a(i) = a(i) - b(j);
        j = j +1;
    end
end
disp("This is the final allocations: ");
disp(t);
disp("The total intials cost is: ");
final_cost_mat = t.*cost;
disp(final_cost_mat);
total_cost = sum(final_cost_mat(:));
disp("Total Transportation Cost:");
disp(total_cost);