% simplex method
var = input("Enter the no of variables: ")
constraints = input("Enter the no of constraints: ")

A = input("Enter the matrix A: ")
b = input("Enter the matrix b: ")
c = input("Enter the matrix c: ")

% we need to add the slack variables
extra_mat = eye(constraints);
A = [A,extra_mat];
disp("This is the updated A: ")
disp(A);
% updated c 
c = [c,zeros(1,constraints)];
disp("This is the updated c: ")
disp(c);

% we need cb now which would all be zero
% handling this will depend on how we make it ie row or coln vector
cb = zeros(1,constraints); 
table = [A,b];
disp("Initial simplex table: ")
disp(table);
zj_minus_cj = cb*table(:,1:end-1)-c;  % remeber we need to multiply the column with the cb
iter = 1; % one done here
while any(zj_minus_cj<0)
    iter = iter + 1;
    disp("No of iteration: ")
    disp(iter);
    if iter> nchoosek(var+constraints,constraints)
        disp("No solution exist. ")
        return;
    end
    [most_neg,entering_col] = min(zj_minus_cj);
    % now ratios
    ratio = zeros(constraints,1);
    for i = 1:constraints
        if table(i,entering_col)>0
            ratio(i) = table(i,end)/table(i,entering_col);
        end
        if table(i,entering_col)<=0
            ratio(i) = inf;
        end
    end
    if all(ratio==inf)
        disp("Unbouned solution.");
        return;
    end
    [value,leaving_col] = min(ratio);
    pivot = table(leaving_col,entering_col);
    % update the leaving row
    table(leaving_col,:) = table(leaving_col,:)/pivot ;
    % update rest of rows
    for i = 1:constraints
        if i ~=leaving_col
            table(i,:) = table(i,:)-table(i,entering_col)*table(leaving_col,:);
        end
    end
    disp("Updated simplex table: ")
    disp(table)
    % update zj-cj
    % but before update the cb
    cb(leaving_col) = c(entering_col);
    zj_minus_cj = cb*table(:,1:end-1) - c;
    disp("Updated zj-cj: ")
    disp(zj_minus_cj);
end

% now we will make the solution
% make whole solution
X = zeros(1,var+constraints);
for i = 1 : size(table,2)-1   % excluding the RHS column
    if sum(table(:,i)==1)==1 && sum(table(:,i)~=0)==1
        % means this column is the basic col
        ind = find(table(:,i)==1);
        X(i) = table(ind,end);
    end
end   
disp("Final solution X:")
disp(X);
solution = cb*table(:,end)

