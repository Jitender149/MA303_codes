% two phase self
var = input("Enter the no of variables: ");
constraints = input("Enter the no of constraints: ");
less_eq = input("Enter the no of less than equal to constraints: ");
eq = input("Enter the no of equal to constraints: ");
great_eq = input("Enter the no of greater than equal to constraints: ");
A = input("Enter the matrix A: ");
b = input("Enter the matrix b: ");
c = input("Enter the matrix c: ");


% Phase I:
extra_mat = [];
extra_mat_c = zeros(1,less_eq+eq+2*great_eq);

for i = 1 : less_eq
    col = zeros(constraints,1);
    col(i) = 1;
    extra_mat = [extra_mat,col];
end
for i = 1 : eq
    col = zeros(constraints,1);
    col(i+less_eq) = 1;
    extra_mat = [extra_mat,col];
end
for i = 1 : great_eq
    col1 = zeros(constraints,1);
    col2 = zeros(constraints,1);
    col1(i+less_eq+eq) = -1;
    col2(i+less_eq+eq) =1;
    extra_mat = [extra_mat,col1];
    extra_mat = [extra_mat,col2];
end
A = [A,extra_mat];
disp("This is the updated A: ");
disp(A);
% c would be just -1 for artificial
c_p1 = zeros(1,var);
for i = 1 : eq
    extra_mat_c(i+less_eq) = -1;
end
for i = 1 : great_eq
    extra_mat_c(i+less_eq+eq+1)=-1;
end
c_p1 = [c_p1,extra_mat_c];
disp("This is the phase 1 c: ");
disp(c_p1);
% we need cb and zj-cj now
cb = zeros(1,constraints);
for i = 1 : eq+great_eq
    cb(i+less_eq) = -1;
end
disp("This is the cb: ");
disp(cb);
zj_minus_cj = cb*A-c_p1;
table = [A, b];
disp("Initial simplex table: ");
disp(table);
disp("This is intial zj-cj: ");
disp(zj_minus_cj);
iter = 1;
max_iter = nchoosek(var+constraints,constraints);
% keeping which var are basic variable
basic_var = zeros(1,constraints);
% initially slack and eq and sur artificial are basic var
for i = 1:less_eq
    basic_var(i) = var + i;
end
for i = 1 : eq
    basic_var(i+less_eq) = var + less_eq + i;
end
for i = 1 : great_eq
    basic_var(i+less_eq+eq) = var + less_eq + eq + i +1 ; % as first surplus and then artificial 
end
disp("These are initial basic variables: ");
disp(basic_var);
while any(zj_minus_cj<0)
    iter = iter+1;
    if iter>max_iter
        disp("No solution");
        return;
    end
    [min_value,entering_col] = min(zj_minus_cj);
    ratio = inf(constraints,1);
    for i = 1 : constraints
        if table(i,entering_col)>0
            ratio(i) = table(i,end)/table(i,entering_col);
        end
    end
    if all(ratio==inf)
        disp("Unbounded solution");
        return;
    end
    [min_ratio,pivot_row] = min(ratio);
    pivot = table(pivot_row,entering_col);
    table(pivot_row,:) = table(pivot_row,:)/pivot;

    for i = 1 : constraints
        if i~=pivot_row
            table(i,:) = table(i,:) - table(i,entering_col)*table(pivot_row,:);
        end
    end
    % update cb and zj-cj
    cb(pivot_row) = c_p1(entering_col);
    basic_var(pivot_row) = entering_col; % directly getting the column here
    disp("These are basic variable: ");
    disp(basic_var);
    zj_minus_cj = cb*table(:,1:end-1)-c_p1;
     disp("This is updated table: ");
    disp(table);
    disp("This is updated cb: ");
    disp(cb);
    disp("This is updated zj-cj: ");
    disp(zj_minus_cj);
end
% now at the end of phase artificial variable can be present 
% if possible do pivoting and if not possible then redundancy 
% here I have assumed that redundancy will not be there
% so we need to do pivoting and remove the artificial variable 
% as per told in class pivoting and redundancy not considered

% phase II:
% here cb changes and c changes and rest table remains same except drop the
% artificial columns
artficial_col = find(c_p1 == -1);
disp("These are the artificial variables: ");
disp(artficial_col);
table(:,artficial_col) = [];
disp("This is the updated table: ");
disp(table);
c_p2 = c;
% and add 0 for slack and surplus var
extra_mat_c_1 = zeros(1,less_eq+great_eq);
c_p2 = [c_p2,extra_mat_c_1];
disp("This is the c for phase 2: ");
disp(c_p2);
% we already have basic var above
cb1 = c_p2(basic_var);
disp("This is cb1 for phase 2: ");
disp(cb1);
zj_minus_cj_1 = cb1*table(:,1:end-1)-c_p2;
iter_1 = 1;
max_iter_1 = max_iter;
while any(zj_minus_cj_1<0)
    iter_1 = iter_1 + 1;
    if iter_1>max_iter
        disp("No solution");
        return;
    end
    [most_neg,entering_col] = min(zj_minus_cj_1);
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
    % update cb1 and zj-cj_1
    cb1(pivot_row) = c_p2(entering_col);
    zj_minus_cj_1 = cb1*table(:,1:end-1)-c_p2;
    disp("This is updated table: ");
    disp(table);
    disp("This is updated cb: ");
    disp(cb1);
    disp("This is updated zj-cj: ");
    disp(zj_minus_cj_1);
end
% now making the solution
X_1 = zeros(1,var+constraints);
for i = 1 : size(table,2)-1  % as RHS at last position
   if sum(table(:,i)==1)==1 && sum(table(:,i)~=0)==1
        basic_var_1 = find(table(:,i)==1);
        X_1(i)=table(basic_var_1,end);
   end
end
disp("This is the final solution: ");
disp(X_1);
soln = cb1*table(:,end)