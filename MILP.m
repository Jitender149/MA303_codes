var = input("Enter the no of variables: ");
constraints = input("Enter the no of constraints: ");
less_eq = input("Enter the no of less than equal to constraints: ");
eq = input("Enter the no of equal to constraints: ");
great_eq = input("Enter the no of greater than equal to constraints: ");
A = input("Enter the matrix A: ");
b = input("Enter the matrix b: ");
c = input("Enter the matrix c: ");
int_vars = input("Enter the varibales that are constrained to be integers ");
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
% now once we got the optimal table from Simplex(Big M) method we need to
% represent the table in the form as per method ie in first row we need
% | z | optimal value | zj-cj | 
% so lets make in that form
% so we need optimal value before
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
% so now we will make the row to append
row1 = zeros(1,1);
row1(1,1) = soln;
row1 = [zj_minus_cj,row1];
disp("This is the row to append as row 1 in the table: ");
disp(row1);
% so get the new table
table = [row1;table];
disp("New table: ");
disp(table);

% Now we will get into the ILPP iterations
while true
    disp("OUTER_LOOP");
    % first we will check if the solution obtained above is already ALL
    % integer then good 
    is_int = true;
    for i = 1 : size(table,1)
        if int_idx(i)==1  % to check only for the constrained varibales
        % last column in my table has the xb values
        % we need to check if these all are integers or not
        temp2 = table(i,end) - round(table(i,end));
        %disp("temp2 it is: ");
        %disp(temp2);
        if abs(temp2) >= 1e-5
            is_int = false;
            break;
        end
        end
    end
    if is_int==true
        disp("Integer soln found.")
        % means this is the optimal soln only
        disp("This is the optimal soln: ");
        disp(table(1,end));
        return;
    end
    % else we go and check the fractional parts of the xb and get the max
    % fractional part row and us that to make the Gomory's cut and append
    % that in the table
    max_frac_part = 0;
    max_frac_part_idx = 0;
    for i = 1 : size(table,1)
        if int_vars(i)==1
            frac_part = table(i,end) - floor(table(i,end));
                if frac_part>max_frac_part && frac_part<1
                    max_frac_part_idx = i;
                    max_frac_part = frac_part;
                end
        end
    end
    disp("This is the max fractional part in the current table: ");
    disp(max_frac_part);
    disp("And this is at row: ");
    disp(max_frac_part_idx);
    % now we make the new constraint
    % we need the non basic variables in the current table (R)
    % and now we need to write these in R1, R2+ and R2-
    % R1 are non basic variables which are constrained to be integers
    % R2+ as those from R-R1 which has yij>0 in the max fraction line
    % and R2- are those from R-R1 which has yij<0  
    % What about yij==0 case ? Does these exist at all?
    % here as mixed integers then we dont guarantee that Z will be integers
    % so not in int_idx
    basic_var_indices = [];
    non_basic_var_indices = [];
    for i = 1 : size(table,2)-1
        temp3 = zeros(1,1);
        temp3 = i;
        if sum(table(1:size(table,1),i)==1)==1 && sum(table(1:size(table,1),i)~=0)==1  % first row has to be left as written above
            % means this is a basic col
            basic_var_indices = [basic_var_indices,temp3];
        else
            non_basic_var_indices = [non_basic_var_indices,temp3];
        end
    end
    disp("These are the non basic indices: ");
    disp(non_basic_var_indices);
    disp("This is int_idx_updated: ");
    disp(int_idx);
    R1 = [];
    R2 = [];
    R2_plus = [];
    R2_minus =[];
    % now we need to make R1, R2+ and R2-
    for i = 1 : size(non_basic_var_indices,2)
        temp4 = zeros(1,1);
        temp4 = non_basic_var_indices(i);
        if int_idx(temp4)==1 % ie non basic var is constrained to be integer
            R1 = [R1,temp4];
        else 
            R2 = [R2,temp4];
        end
    end
    disp("R2:");
    disp(R2);
    for i = 1 : size(R2,2)
        temp5 = zeros(1,1);
        temp5 = non_basic_var_indices(i);
        disp("This is element: ");
        disp(table(max_frac_part_idx,temp5));
        if table(max_frac_part_idx,temp5)>0
            R2_plus = [R2_plus,temp5];
        else
            R2_minus = [R2_minus,temp5];
        end
    end
    % now we have R1, R2_plus, R2_minus
    disp("R1 ");
    disp(R1);
    disp("R2_plus ");
    disp(R2_plus);
    disp("R2_minus ");
    disp(R2_minus);
    new_row = zeros(1,size(table,2)+1);  % new cut as 1 more var thus table col also increases by 1
    disp("Debugging: ");
    disp(table(max_frac_part_idx,end));
    disp(floor(table(max_frac_part_idx,end)))
    new_row(1,end) = -1*(table(max_frac_part_idx,end)-floor(table(max_frac_part_idx,end)));
    new_row(1,end-1) = 1; % new variable ke corresponding ka 1
    % first corresponding to R1
    for i = 1: size(R1,2)
        idx = non_basic_var_indices(1,i);
        new_row(1,idx) = -1*(table(max_frac_part_idx,idx) - floor(table(max_frac_part_idx,idx)));
    end

    for i = 1 : size(R2_plus,2)
        idx  = non_basic_var_indices(1,R2_plus(i));
        new_row(1,idx) = -1*(table(max_frac_part_idx,idx));
    end

    for i = 1 : size(R2_minus)
        idx = non_basic_var_indices(1,R2_minus(i));
        fractionnnn = table(max_frac_part_idx,idx) - floor(table(max_frac_part_idx,idx));
        new_row(1,idx) = fractionnnn*table(max_frac_part_idx,idx)/(1-fractionnnn);
    end

    retarded_col = zeros(size(table,1),1);
    table = [table(:,1:end-1),retarded_col,table(:,end)];
    table = [table;new_row];
    disp("This is table after appending the Gomory's cut: ");
    disp(table);
    % and c would also change with this, so would cb 
    c = [c,zeros(1,1)];
    disp("The updated c: ");
    disp(c);
    cb = [cb,zeros(1,1)];
    disp("The updated cb: ");
    disp(cb);
    % Now we would apply dual simplex method to get next feasible table
        % Now we would apply dual simplex method to get next feasible table
    iter = 1;
    max_iter = 100; % safety limit
    while true
        disp("INNER_LOOP");
        iter = iter + 1;
        if iter > max_iter
            disp("Max iterations reached in dual simplex");
            return;
        end
        % Step 1: Identify row with most negative b_i
        b_column = table(2:end,end); % skip row 1 (zj-cj)
        [min_bi, leaving_row_relative] = min(b_column);
        if min_bi >= 0
            disp("Feasible solution obtained from dual simplex.");
            break;  % stop dual simplex
        end
        leaving_row = leaving_row_relative + 1; % offset due to first row

        % Step 2: Identify pivot column (min ratio of |zj-cj / a_ij| where a_ij < 0)
        ratios = inf(1, size(table,2)-1); % last column is b, not needed
        for j = 1:size(table,2)-1
            if table(leaving_row, j) < 0
                ratios(j) = abs(table(1,j) / table(leaving_row, j));
            end
        end

        [min_ratio, entering_col] = min(ratios);
        if min_ratio == inf
            disp("Dual simplex failed: No valid pivot column (unbounded).");
            return;
        end

        % Step 3: Pivot operation
        pivot = table(leaving_row, entering_col);
        table(leaving_row,:) = table(leaving_row,:) / pivot;

        for i = 1:size(table,1)
            if i ~= leaving_row
                table(i,:) = table(i,:) - table(i,entering_col) * table(leaving_row,:);
            end
        end

        % Step 4: Update cb and zj-cj
        cb(leaving_row-1) = c(entering_col); % offset since table's row 1 is zj-cj
        table(1,1:end-1) = cb * table(2:end,1:end-1) - c;  % recompute zj-cj
        table(1,end) = cb * table(2:end,end);  % recompute z
        disp("Dual simplex step completed. Updated table:");
        disp(table);
    end
end