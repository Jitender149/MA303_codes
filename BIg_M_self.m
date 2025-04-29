% big m method
var = input("Enter the no of variables: ");
constraints = input("Enter the no of constraints: ");
less_eq = input("Enter the no of less than equal to constraints: ");
eq = input("Enter the no of equal to constraints: ");
great_eq = input("Enter the no of greater than equal to constraints: ");
A = input("Enter the matrix A: ");
b = input("Enter the matrix b: ");
c = input("Enter the matrix c: ");

extra_mat = [];
extra_mat_c = zeros(1,less_eq+eq+2*great_eq);
% for less_eq --> slack
for i = 1:less_eq
    col = zeros(constraints,1);
    col(i) = 1;
    extra_mat = [extra_mat,col];
end

% for eq --> artificial 
for i = 1:eq
    col = zeros(constraints,1);
    col(i+less_eq) = 1;
    extra_mat = [extra_mat,col];
end

% for great_eq --> surplus and artificial
for i = 1:great_eq
    col1 = zeros(constraints,1);
    col2 = zeros(constraints,1);
    col1(i+less_eq+eq) = -1;
    col2(i+less_eq+eq) = 1;
    extra_mat = [extra_mat,col1];
    extra_mat = [extra_mat,col2];
end
% updating extra_mat_c
for i = less_eq+1:less_eq+eq
    extra_mat_c(i) = -1000;
end
for i = less_eq+eq+1:less_eq+eq+great_eq
    extra_mat_c(i+1) = -1000;
end
A = [A, extra_mat];
disp("This is updated A: ")
disp(A);
c = [c,extra_mat_c];
disp("This is updated c: ")
disp(c);

cb = zeros(1,constraints);
for i = less_eq+1:less_eq+eq+great_eq
    cb(i) = -1000;
end
disp("This is cb: ");
disp(cb);
% now continue the normal simplex
table = [A, b];
disp("This is initial simplex table: ");
disp(table);
zj_minus_cj = cb*A-c;
disp("This is initial zj-cj: ")
disp(zj_minus_cj);
iter = 1;
max_iter = nchoosek(var+constraints,constraints);
while any(zj_minus_cj<0)
    iter = iter + 1;
    if iter>max_iter
        disp("No solution");
        return;
    end
    [most_neg,entering_col] = min(zj_minus_cj);
    ratio = zeros(constraints,1);
    for i = 1:constraints
        if table(i,entering_col)>0
            ratio(i) = table(i,end)/table(i,entering_col);
        else
            ratio(i) = inf;
        end
    end
    if all(ratio==inf)
        disp("Unbounded solution");
        return;
    end
    [min_ratio,pivot_row] = min(ratio);
    pivot = table(pivot_row,entering_col);
    table(pivot_row,:) = table(pivot_row,:)/pivot;

    for i = 1:constraints
        if i~=pivot_row
            table(i,:) = table(i,:) - table(i,entering_col)*table(pivot_row,:);
        end
    end
    % update cb and zj-cj
    cb(pivot_row) = c(entering_col);
    zj_minus_cj = cb*table(:,1:end-1)-c;
    disp("This is updated table: ");
    disp(table);
    disp("This is updated cb: ");
    disp(cb);
    disp("This is updated zj-cj: ");
    disp(zj_minus_cj);
end
X = zeros(1,var+constraints);
for i = 1: size(table,2)-1
    if sum(table(:,i)==1)==1 && sum(table(:,i)~=0)==1
       ind = find(table(:,i)==1); % means jth col is the basic var and its values is in the end of this ind th row
       X(i) = table(ind,end);
    end
end
disp("Final soln:");
disp(X(1:var)); % only the original variables
% optimal answer
soln = cb*table(:,end)