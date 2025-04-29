% Algebraic method 
% taking the inputs
A = input("Enter the matrix A: ");
b = input("Enter the matrix b: ");
c = input("Enter the matrix c: ");

[m,n] = size(A)   % remember this order
disp(n);
disp(m);
% now we need to produce all the combinations 
combinations = nchoosek(1:n,m); 
disp("These are all the combinations: ");
disp(combinations);

feasSol = []; % initialise empty to store all the feasible solutions
% now for each combination make the submatrix and get the solution
for i = 1:size(combinations,1)
    selecCol = combinations(i,:) 
    subMat = A(:,selecCol);
    
    if det(subMat)~=0
    %subMat = inv(subMat)
    soln = subMat\b;
    disp("This is the submatrix: ");
    disp(subMat);
    % now check if this basic feasible solutions
        if all(soln>=0)
            % then make full solution and then append to the feasMat
            x = zeros(n,1);
            for j = 1:m
                x(selecCol(j)) = soln(j);
            end
            disp("This is the solution");
            disp(x);
            % now append to the feasSol matrix
            feasSol = [feasSol,x];
        end
    end
end
disp("This is the feasSol matrix: ");
disp(feasSol);
% to get the max answer we multiply by c and get max from this
maxMat = c*feasSol
[maxi,index] = max(maxMat);
disp("This is the max soln: ")
disp(maxi)
disp("This is the index of solution in the feas matrix: ")
disp(index);
% done and dusted
