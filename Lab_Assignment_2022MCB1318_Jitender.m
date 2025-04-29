% 2022MCB1318
% Jitender Jangra
% Lab Assignment (Assignment Problem)

% Taking input
% Enter large value say 1000 (ofc this is contextual but works here)
% As the question given is balanced only ie square matrix so I didnt handle
% the case explicilty
% though we can add 0 column/row while inputting the matrix and the code
% works for the general assignment problem
cost = input("Enter the cost matrix: ");
flag = input("Enter 1 if maximisation problem else 0: ");
[m,n] = size(cost);
A = cost;  % ofc we dont want to loose original matrix
if flag==1
    % this is maximisation, then we need to select the max element and
    % subtract all elements from this except the Large Value cells as those
    % are BLOCKED places
    sorted_elements = sort(A(:));
    max_element = sorted_elements(n*m-2);  % though this is second max actaully
    disp(max_element);
    % but large values is not in the play so we take second max here
    for i = 1 : m
        for j = 1 : n
            if cost(i,j)~=1000
                A(i,j) = max_element - A(i,j);
            end
        end
    end
end
disp("The input cost matrix is: ");
disp(cost);
disp("This is the cost matrix we will be working with: ");
disp(A);
% now we will start with the method taught in the class...tbh idk the name
% of the method
%Calculating Row Minimum Matrix
for i = 1:n
    A(i,:) = A(i,:) - min(A(i,:));
end
disp("Matrix after row operation: ");
disp(A);
%Calculating Column Minimum Matrix
for i = 1:n
    A(:,i) = A(:,i) - min(A(:,i));
end
disp("Matrix after col operation: ");
disp(A);
% Now we will do the iterations to cover the matrix with line 
% by using method taught in class
while(1)
    disp("Hello! Master JJ here.")
    taken_zero = zeros(n,n); % this wil track the assingned cells (ofc in this method 0 values cells will be assigned)
    vis_row = zeros(n,1); % this will track the marked columns
    vis_col = zeros(n,1); % this will track the marked columns
    % Getting the Initial Assignment
    while(1)
        cnt = 0;
        cnt_row = zeros(n,1);
        cnt_col = zeros(n,1);
        for i = 1:n
            for j = 1:n
                if(A(i,j)==0 && vis_row(i)==0 && vis_col(j)==0)
                    cnt_row(i) = cnt_row(i) + 1;
                    cnt_col(j) = cnt_col(j) + 1;
                end
            end
        end
        row = -1;
        col = -1;
        for i = 1:n
            for j = 1:n
                if(A(i,j)==0 && vis_row(i)==0 && vis_col(j)==0)
                    if(cnt_row(i)==1 || cnt_col(j)==1)
                        row = i;
                        col = j;
                    else 
                        if(row == -1)
                            row = i;
                            col = j;
                        end
                    end
                end
            end
        end
        if(row == -1)
            break
        end
        taken_zero(row,col) = 1;
        vis_row(row) = 1;
        vis_col(col) = 1;
        cnt = cnt + 1;
    end
    %Marking all not-taken zeros
    for i = 1:n
        for j = 1:n
            if(A(i,j)==0 && taken_zero(i,j)==0)
                taken_zero(i,j)=2;
            end
        end
    end
    %now we will tick the non-visited rows
    tick_row = zeros(n,1);
    tick_col = zeros(n,1);
    for i = 1:n
        if(vis_row(i)==0)
            tick_row(i) = 1;
        end
    end
    %now we will tick the columns and rows according to the rules
    % tick the non assigned rows
    % tick the columns which have corssed zero in the marked rows
    % tick the rows (not already marked ofc) which have assignmnts in
    % marked columns
    % this while loop below handles last two steps
    while(1)
        cnt = 0;
        for i = 1:n
            if(tick_row(i) == 1)
                for j = 1:n
                    if(taken_zero(i,j)==2 && tick_col(j)==0)
                        cnt = cnt + 1;
                        tick_col(j) = 1;
                    end
                end
            end
        end
        for j = 1:n
            if(tick_col(j) == 1)
                for i = 1:n
                    if(taken_zero(i,j)==1 && tick_row(i)==0)
                        cnt = cnt + 1;
                        tick_row(i) = 1;
                    end
                end
            end
        end
        if(cnt == 0)
            break
        end
    end
    %now we will Draw and Count Lines
    cnt_lines = 0;
    for i = 1:n
        if(tick_row(i)==0)
            cnt_lines = cnt_lines + 1;
        end
        if(tick_col(i)==1)
            cnt_lines = cnt_lines + 1;
        end
    end
    if(cnt_lines == n) % ofc this has be to satisfied
        break
    end
    %Now perform the operation for minimum element to make more zeros
    % ie we will subtract this min element from 
    % all the non covered elements
    % and add to the elemtns which lie on intersection of two lines
    % rest elements remains unchanged
    min_ele = 1000000000;
    for i = 1:n
        for j = 1:n
            if(tick_row(i) == 1 && tick_col(j) == 0)
                min_ele = min(min_ele,A(i,j));
            end
        end
    end
    for i = 1:n
        for j = 1:n
            if(tick_row(i) == 1 && tick_col(j) == 0)
                A(i,j) = A(i,j) - min_ele;
            end
            if(tick_row(i)==0 && tick_col(j) == 1)
                A(i,j) = A(i,j) + min_ele;
            end
        end
    end
end
taken_zero = zeros(n,n);
vis_row = zeros(n,1);
vis_col = zeros(n,1);
cnt = 0;
%Now, Calculating the final Assignment
while(cnt<n)
    cnt_row = zeros(n,1);
    cnt_col = zeros(n,1);
    for i = 1:n
        for j = 1:n
            if(A(i,j)==0 && vis_row(i)==0 && vis_col(j)==0)
                cnt_row(i) = cnt_row(i) + 1;
                cnt_col(j) = cnt_col(j) + 1;
            end
        end
    end
    row = -1;
    col = -1;
    for i = 1:n
        for j = 1:n
            if(A(i,j)==0 && vis_row(i)==0 && vis_col(j)==0)
                if(cnt_row(i)==1 || cnt_col(j)==1)
                    row = i;
                    col = j;
                else 
                    if(row == -1)
                        row = i;
                        col = j;
                    end
                end
            end
        end
    end
    if(row == -1)
        break
    end
    taken_zero(row,col) = 1;
    vis_row(row) = 1;
    vis_col(col) = 1;
    cnt = cnt + 1;
end
% Now we have our optimal assingment with us and we will calculate the
% cost/value using the orignal matrix and assignments obtained from above
opt_val = 0;
for i = 1:n
    for j = 1:n
        opt_val = opt_val + (taken_zero(i,j)) * cost(i,j);
    end
end
disp("The optimal value is: ");
disp(opt_val);
disp("Assignment is: ")
disp(taken_zero);
% now looking from this we see that: 
disp("Optimal Assignment:");
for i = 1:n
    for j = 1:n
        if taken_zero(i,j) == 1
            fprintf("Pilot %d → Flight No %d (Rating = %d)\n", i, j, cost(i,j));
        end
    end
end