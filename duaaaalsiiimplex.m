% dual simplex method
% solving assuming that we are solving maximisation problem as this is how
% taught in class
% if minimisation problem then just take -ve of final answer
% first we need to take input
var = input("Enter the no of variables: ");
constraints = input("Enter the no of constraints: ");
less_eq = input("Enter the no of less than equal to constraints: ");
great_eq = input("Enter the no of greater than equal to constraints: ");
eq = input("Enter the no of equal to constraints: ");
% give in order in which writing type of constraint
A = input("Enter the matrix A: ");
b = input("Enter the matrix b: ");
c = input("Enter the matrix c: ");

% for less than equal to constraints we simply add slack variable
% for equal to constraints we break it into less than and greater than
% for greater than type constraints we simply add surplus variable

extra_mat = [];
extra_mat_c = zeros(1,less_eq+great_eq+2*eq);
A_new = [];
temp = [];

for i = 1:less_eq
    col = zeros(constraints+eq,1);
    col(i) = 1;
    extra_mat = [extra_mat,col];
    A_new = [A_new;A(i,:)];
    temp = [temp,col];
end

disp("This is A_new: ");
disp(A_new);
disp("This is temp: ");
disp(temp);

for i = 1:great_eq
    col1 = zeros(constraints+eq,1);
    col1(i+less_eq) = -1;
    extra_mat = [extra_mat,col1];
    A_new = [A_new;-1*(A(i+less_eq,:))];
    temp = [temp,-1*(col1)];
end

disp("This is A_new: ");
disp(A_new);
disp("This is temp: ");
disp(temp);

for i = 1:eq
    col1 = zeros(constraints+eq,1);
    col1(i+less_eq+great_eq) = 1;
    extra_mat = [extra_mat,col1];
    A_new = [A_new;A(i+less_eq+great_eq)];
    temp = [temp,col1];
end

for i = 1:eq
    col2 = zeros(constraints+eq,1);
    col2(i+less_eq+great_eq+eq) = -1;
    extra_mat = [extra_mat,col2];
    A_new =[A_new;-1*(A(i+less_eq+great_eq))];
    temp = [temp,-1*(col2)];
end

disp("This is temp");
disp(temp);
A_new = [A_new,temp];
A = A_new;

for i = 1:less_eq+2*eq+great_eq
    extra_mat_c(i) = 0;
end

disp("This is updated A: ");
disp(A);
c = [c,extra_mat_c];
disp("This is updated c: ");
disp(c);
cb = extra_mat_c;
disp("This is cb: ");
disp(cb);

b_new = zeros(constraints+eq,1);
for i = 1 : less_eq
    b_new(i) = b(i);
end
for i = 1 : great_eq
    b_new(less_eq+i) = -1*b(less_eq+i);
end
for i = 1 : eq
    b_new(i+less_eq+great_eq) = b(i + less_eq+great_eq);
end
for i = 1 : eq
    b_new(i+less_eq+great_eq+eq) = -1*b(i + less_eq+ great_eq);
end

disp("This is b_new: ");
disp(b_new);
table = [A,b_new];
disp("Initial simplex table: ")
disp(table);
zj_minus_cj = cb*table(:,1:end-1)-c;
disp("Initial zj-cj: ");
disp(zj_minus_cj);
iter = 1;
while any(zj_minus_cj<0)
    iter = iter + 1;
    disp("No of iteration: ");
    disp(iter);
    if iter> nchoosek(var+constraints,constraints)
        disp("Bhagwan jaane iska most prolly non optimality ")
        return;
    end
    [most_neg,entering_col] = min(zj_minus_cj);
    ratio = zeros(constraints+eq,1);
    for i = 1:constraints+eq
        if table(i,entering_col)>0
            ratio(i) = table(i,end)/table(i,entering_col);
        end
        if table(i,entering_col)<=0
            ratio(i) = inf;
        end
    end
    if all(ratio==inf)
        disp("Unbouned solution....bhagwaan jaane iska bhi bhai.");
        return;
    end
    [value,leaving_col] = min(ratio);
    pivot = table(leaving_col,entering_col);
    table(leaving_col,:) = table(leaving_col,:)/pivot;
    for i = 1:constraints+eq
        if i ~=leaving_col
            table(i,:) = table(i,:)-table(i,entering_col)*table(leaving_col,:);
        end
    end
    disp("Updated simplex table: ")
    disp(table)
    cb(leaving_col) = c(entering_col);
    zj_minus_cj =cb*table(:,1:end-1)-c;
    disp("Updated zj-cj: ")
    disp(zj_minus_cj);
end

while any(table(:,end)<0)
    disp("Dual Simplex Iteration");
    [most_neg_xb, leaving_row] = min(table(:,end));
    ratio = -inf(1,size(table,2)-1); % ie all -inf
    cntr = 0;
    for i = 1 : size(table,2)-1
        if table(leaving_row,i)<0
        ratio(i) = zj_minus_cj(i)/table(leaving_row,i);
        else
            cntr = cntr + 1; 
        end
    end
    if cntr==size(table,2)-1
        disp("Infeasible solution")
        return
    end
    disp("These are the ratios: ");
    disp(ratio);
    [max_ratio, entering_col] = max(ratio);
    pivot = table(leaving_row, entering_col)
    table(leaving_row,:) = table(leaving_row,:) / pivot;
    for i = 1:constraints+eq
        if i ~= leaving_row
            table(i,:) = table(i,:) - table(i,entering_col) * table(leaving_row,:);
        end
    end
    disp("Updated simplex table after dual simplex iteration: ");
    disp(table);
    cb(leaving_row) = c(entering_col);
    disp("The updated cb: ");
    disp(cb);
    zj_minus_cj = cb*table(:,1:end-1) - c;
    disp("Updated zj-cj after dual simplex iteration: ");
    disp(zj_minus_cj);
end

disp("Optimal solution reached!");
disp("Final Simplex Table: ");
disp(table);
ans = cb*table(:,end);
disp("Optimal solution: ");
disp(ans);
% handle if the minimisation problem i just consider -ve of the above
% optimal wrt maximisation
%ratio = abs(zj_minus_cj ./ table(leaving_row,1:end-1)) % only if 
%ratio(table(leaving_row,1:end-1) >= 0) = inf