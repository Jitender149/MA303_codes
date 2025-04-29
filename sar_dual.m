clc
nv = input("Enter the number of variables: ");
nc = input("Enter the number of constraints: ");
L = input("Enter the no. of less than equal to constraints: ");
E = input("Enter the no. of equal to constraints: ");
G = input("Enter the no. of greater than equal to constraints: ");
A = input("Enter the matrix A: ");
b = input("Enter the constant matrix (RHS): ");
c = input("Enter the coefficients of the objective function: ");
C = c;
bv = zeros(0);
for i = 1:L
    t = zeros([nc,1]);
    t(i,1) = 1;
    A = [A t];
    C(nv+i,1) = 0;
    bv(i,1) = nv+i;
end
for i = 1:G
    t = zeros([nc,1]);
    [x,y] = size(A);
    A(i+L,1:y) = -1*A(i+L,1:y);
    t(i+L,1) = 1;
    b(i+L,1) = -1*b(i+L,1);
    A = [A t];
    C(nv+L+i,1) = 0;
    bv(i+L,1) = nv+L+i;
end
A
A1 = zeros(0);
b1 = zeros(0);
for i = 1:L+G
    [x,y] = size(A);
    A1 = [A1' A(i,1:y)']';
    b1 = [b1 b(i,1)];
end
for i = 1:E
    [x,y] = size(A);
    A1 = [A1' A(i+L+G,1:y)']';
    A1 = [A1' A(i+L+G,1:y)']';
     b1 = [b1 b(i+L+G,1)];
     b1 = [b1 b(i+L+G,1)];
end
A  = A1;

b= b1';
 [x,y] = size(A);
 nc = x;
for i = 1:E
    t = zeros([nc,1]);
    t(i+L+G,1) = 1;
    A = [A t];
    C(nv+i+L+G,1) = 0;
    bv(i+L+G,1) = nv+L+i+G;
end
for i = 1:E
    t = zeros([nc,1]);
    [x,y] = size(A);
    A(i+L+G+E,1:y) = -1*A(i+L+G+E,1:y);
    t(i+L+G+E,1) = 1;
    b(i+L+G+E,1) = -1*b(i+L+G+E,1);
    A = [A t];
    C(nv+L+i+G+E,1) = 0;
    bv(i+L+G+E,1) = nv+L+i+G+E;
end
A
b
C
bv
[nc,n]= size(A);
k = nchoosek(n,nc);
y =0;
co = zeros(0);
for i = 1:L+2*E+G
    co(i,1) = C(bv(i,1),1);
end
table = [bv co b A];
no = zeros(0);
no(1,1) =0;
no(2,1) = 0;
no(3,1) = 0;
for i = 1:n
    no(i+3,1) = i;
end
table = [no table']';
table

while(y<k)
    for i = 1:n
        t = 0;
        for j = 1:nc
            t= t+table(j+1,2)*table(j+1,i+3);
        end
        z_c(1,i) = t - C(i,1);
    end
    z_c
    t =0;
    for i =1:n
        if(z_c(1,i)<0)
            t=t+1;
        end
    end
    if(t==0)
         disp("piovting done");
        break;
    end
    t = 0;
    for i = 1:n
    if(t>z_c(1,i))
        pc = i+3;
        t =z_c(1,i);
    end
    end
    t =0;
    for i = 1:nc
        if(table(i+1,pc)>0)t =1;
        end
    end
    if(t==0)
        disp("infeasible solution");
        return;
    end
t = 500000;
for i=1:nc
    if(table(i+1,pc)>0)
        u = table(i+1,3)/table(i+1,pc);
        if(u<=t)
            pr= i+1;
            t = u;
        end
    end
end
pc
pr
ntable = table;
ntable(pr,1) = table(1,pc);
ntable(pr,2) = C(pc-3,1);

for i=1:n+1
    ntable(pr,i+2) = ntable(pr,i+2)/table(pr,pc);
end
for j=2:nc+1
    if(j==pr)
        continue;
    end
    for i = 1:n+1
        ntable(j,i+2) = ntable(j,i+2)-ntable(pr,i+2)*table(j,pc);
    end
end
table  = ntable;
table
y = y+1;
end
table
k = nchoosek(n,nc);
y =0;
while(y<k)
    for i = 1:n
        t = 0;
        for j = 1:nc
            t= t+table(j+1,2)*table(j+1,i+3);
        end
        z_c(1,i) = t - C(i,1);
    end
    z_c
    t =0;
    for i =1:nc
        if(table(i+1,3)<0)
            t=t+1;
        end
    end
    if(t==0)
        disp("solution finded");
        for i = 1:n
            sol(i,1) =0;
        end
        for i =1:nc
            sol(table(i+1,1),1) = table(i+1,3);
        end
        sol
        disp('value is : ')
        (C')*sol
        return;
    end
    t = 0;
    for i = 1:nc
    if(t>table(i+1,3))
        pr = i+1;
        t =table(i+1,3);
    end
    end
    t =0;
    for i = 1:n
        if(table(pr,i+3)<0)t =1;
        end
    end
    if(t==0)
        disp("infeasible solution");
        return;
    end
t = -500000;
for i=1:n
    if(table(pr,i+3)<0)
        u = z_c(1,i)/table(pr,i+3);
        if(u>=t)
            pc= i+3;
            t = u;
        end
    end
end
pc
pr
ntable = table;
ntable(pr,1) = table(1,pc);
ntable(pr,2) = C(pc-3,1);

for i=1:n+1
    ntable(pr,i+2) = ntable(pr,i+2)/table(pr,pc);
end
for j=2:nc+1
    if(j==pr)
        continue;
    end
    for i = 1:n+1
        ntable(j,i+2) = ntable(j,i+2)-ntable(pr,i+2)*table(j,pc);
    end
end
table  = ntable;
table
y = y+1;
end