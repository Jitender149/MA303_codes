clc;
clear;
% Example Input (Replace with your ILP)
% Maximize z = 3x + 2y
% Subject to:
%   2x + y ≤ 4
%   x + 2y ≤ 5
%   x, y ≥ 0 and integers

A = [7 16; 3 -2];
b = [52; 9];
C = [15; 32];
[nc, n] = size(A);


table = zeros(nc+1, n+nc+3);
table(1,4:3+n) = 1:n;
table(1,1:3) = [0 0 0];

for i = 1:nc
    table(i+1,1) = n+i;       
    table(i+1,2) = 0;         
    table(i+1,3) = b(i);    
    table(i+1,4:3+n) = A(i,:);
    table(i+1,3+i+n) = 1;   
end
C = [C; zeros(nc,1)];
y = 0;
k = 50;
table
n =4;
while(y < k)
    for i = 1:n
        t = 0;
        for j = 1:nc
            t = t + table(j+1,2) * table(j+1,i+3);
        end
        z_c(1,i) = t - C(i,1);
    end
    z_c

    t = 0;
    for i = 1:n
        if(z_c(1,i) < 0)
            t = t + 1;
        end
    end

    if(t == 0)
        disp("Initial LP Relaxation done (Simplex Pivoting Finished)");
        break;
    end

    t = 0;
    for i = 1:n
        if(t > z_c(1,i))
            pc = i+3;
            t = z_c(1,i);
        end
    end

    t = 0;
    for i = 1:nc
        if(table(i+1,pc) > 0)
            t = 1;
        end
    end

    if(t == 0)
        disp("Unbounded solution");
        break;
    end

    t = 500000;
    for i = 1:nc
        if(table(i+1,pc) > 0)
            u = table(i+1,3) / table(i+1,pc);
            if(u <= t)
                pr = i+1;
                t = u;
            end
        end
    end

    pc
    pr

    ntable = table;
    ntable(pr,1) = table(1,pc);
    ntable(pr,2) = C(pc-3,1);

    for i = 1:n+1
        ntable(pr,i+2) = ntable(pr,i+2) / table(pr,pc);
    end

    for j = 2:nc+1
        if(j == pr)
            continue;
        end
        for i = 1:n+1
            ntable(j,i+2) = ntable(j,i+2) - ntable(pr,i+2) * table(j,pc);
        end
    end

    table = ntable;
    table
    y = y + 1;
end
%ILPP
while true
    
    is_integer = true;
    for i = 1:nc
        val = table(i+1,3);
        if abs(val - round(val)) > 1e-5
            is_integer = false;
            break;
        end
    end

    if is_integer
        disp("Integer solution found:");
        sol =0;
        for i = 1:nc
            
                sol = sol + table(i+1,2)*table(i+1,3);
           
        end
       
        disp("Objective Value:");
        disp(sol);
        return;
    end

    yy =0;
    for i = 1:nc
        if abs(table(i+1,3) - floor(table(i+1,3))) >yy
            frac_row = i+1;
            yy = abs(table(i+1,3) - floor(table(i+1,3)));
        end
    end

    new_row = zeros(1, size(table,2));
    new_row(1) = size(table,2)-2;
   % C = [C' 0]';
    new_row(2) = 0;
    new_row(3) = -mod(table(frac_row,3),1);  % RHS

    for i = 4:size(table,2)
        new_row(i) = -mod(table(frac_row,i),1);
    end

    table = [table; new_row]; % add cut
    table = [table, [size(table,2)-2; zeros(size(table,1)-2,1); 1]]; 
    new_row = [new_row, 1]; 
 

    nc = nc + 1;
    C = [C; 0];

    k = 50;
    y = 0;
    table
    [nc,n] = size(table);
    n = n-3;
    nc = nc-1;
    while(y < k)
        for i = 1:n
            t = 0;
            for j = 1:nc
                t = t + table(j+1,2) * table(j+1,i+3);
            end
            z_c(1,i) = t - C(i,1);
        end
        z_c

        t = 0;
        for i = 1:nc
            if(table(i+1,3) < 0)
                t = t + 1;
            end
        end
        if(t == 0)
            disp("Dual Simplex: Feasible again");
            break;
        end

        t = 0;
        for i = 1:nc
            if(t > table(i+1,3))
                pr = i+1;
                t = table(i+1,3);
            end
        end

        t = 0;
        for i = 1:n
            if(table(pr,i+3) < 0)
                t = 1;
            end
        end
        if(t == 0)
            disp("Infeasible solution");
            return;
        end

        t = -500000;
        for i = 1:n
            if(table(pr,i+3) < 0)
                u = z_c(1,i) / table(pr,i+3);
                if(u >= t)
                    pc = i+3;
                    t = u;
                end
            end
        end

        pc
        pr

        ntable = table;
        ntable(pr,1) = table(1,pc);
        ntable(pr,2) = C(pc-3,1);

        for i = 1:n+1
            ntable(pr,i+2) = ntable(pr,i+2) / table(pr,pc);
        end

        for j = 2:nc+1
            if(j == pr)
                continue;
            end
            for i = 1:n+1
                ntable(j,i+2) = ntable(j,i+2) - ntable(pr,i+2) * table(j,pc);
            end
        end

        table = ntable;
        table
        y = y + 1;
    end
end