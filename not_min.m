A1 = input('Enter matrix for <= type constarint: ');
A2 = input('Enter matrix for = type constraints: ');
A3 = input('Enter matrix for >= type constraint: ');
b = input('Enter the constant term vector for the constraints: ');
c = input('Enter coefficients of objective function: ');
x = input('Enter number of varaibles:');
u = x;

[m1,n1] = size(A1); [m2,n2] = size(A2); [m3,n3] = size(A3);
A = zeros(m1+m2+m3,x+m1+m2+2*m3);
if m1>0
    A(1:m1,1:n1) = A1;
end
if m2>0
    A(m1+1:m1+m2,1:n2) = A2;
end
if m3>0
    A(m1+m2+1:m1+m2+m3,1:n3) = A3;
end

if m1>0
    A(1:m1,x+1:x+m1) = eye(m1);
end
if m2>0
    A(m1+1:m1+m2,x+m1+1:x+m1+m2) = eye(m2);
end
if m3>0
    A(m1+m2+1:m1+m2+m3,x+m1+m2+1:x+m1+m2+m3) = eye(m3);
    A(m1+m2+1:m1+m2+m3,x+m1+m2+m3+1:x+m1+m2+2*m3) = -eye(m3);
end
C = zeros(x+m1+m2+2*m3,1);
C(x+m1+1:x+m1+m2+m3) = -1;
[m,n] = size(A);
S = zeros(m,n+2);
S(m1+1:m1+m2+m3,1) = -1;
S(:,2) = b;
S(:,3:n+2) = A(:,1:n);
disp('Initial Phase 1 Table:');
disp(S);
while true
    ent = 3; p = 0;
    for i = 1:m
        p = p + S(i,1)*S(i,3);
    end
    p = p - C(1);
    for j = 4:n+2
        p1 = 0;
        for i = 1:m
            p1 = p1 + S(i,1)*S(i,j);
        end
        p1 = p1-C(j-2);
        if p1<p
            p = p1;
            ent = j;
        end
    end
    if p >= 0
        disp('Phase 1 Ends');
        break;
    end
    lev = 1; rat = 1000000000;
    for i = 1:m
        if(S(i,ent) > 0 && rat > S(i,2)/S(i,ent))
            rat = S(i,2)/S(i,ent);
            lev = i;
        end
    end
    piv = S(lev,ent);
    S(lev,2:n+2) = S(lev,2:n+2)/piv;
    S(lev,1) = C(ent-2);
    for i = 1:m
        if i ~= lev
            t = S(i,ent);
            S(i,2:n+2) = S(i,2:n+2) - t*S(lev,2:n+2);
        end
    end
end
disp(S)
for i = 1:m
    if S(i,1) < 0
        red = [];
        if S(i,2) ~= 0
            error('Infeasible Solution');
        else
            cnt = 0;
            for j = 3:n+2
                if S(i,j) == 0
                    cnt = cnt+1;
                end
            end
            if cnt == n-1
                red = [red,i];
            else
                for j = 3:n+2
                    ele = S(i,j);
                    if ele ~= 0 && ele ~= 1 
                        S(i,2:n+2) = S(i,2:n+2)/ele;
                        S(i,1) = C(j-2);
                        for k = 1:m
                            if k ~= i
                                t = S(k,j);
                                S(k,2:n+2) = S(k,2:n+2) - t*S(i,2:n+2);
                            end
                        end
                        break;
                    elseif S(i,j) == 1
                        cnt = 0;
                        for h = 1:m
                            if S(h,j) == 0
                                cnt = cnt+1;
                            end
                        end
                        if cnt ~= m-1
                            S(i,2:n+2) = S(i,2:n+2)/ele;
                            S(i,1) = C(j-2);
                            for k = 1:m
                                if k ~= i
                                    t = S(k,j);
                                    S(k,2:n+2) = S(k,2:n+2) - t*S(i,2:n+2);
                                end
                            end
                            break;
                        end
                    end
                end
            end
        end
    end
end
disp(S)
S(:,x+m1+3:x+m1+m2+m3+2) = [];
disp('First Table for Phase-2:');
W = eye(m,m);
n = n-m2-m3;
C = zeros(x+m1+m3);
C(1:x) = c;
for i = 1:m
    for j = 3:n+2
        if isequal(W(:,i),S(:,j))
            S(i,1) = C(j-2);
            break;
        end
    end
end
disp(S)
while true
    p = 0; enter = 3;
    for i = 1:m
        p = p +S(i,1)*S(i,3);
    end
    p = p - C(1);
    for j =4:n+2
        p1 = 0;
        for i = 1:m
            p1 = p1 + S(i,1)*S(i,j);
        end
        p1 = p1-C(j-2);
        if(p1 < p)
            p = p1;
            enter = j;
        end
    end
    if p >= 0
        disp('Optimal Found');
        ansi = 0;
        for i =1:m
            ansi = ansi + S(i,1)*S(i,2);
        end
        disp(ansi);
        break;
    end
    cnt = 0;
    for i = 1:m
        if(S(i,enter) < 0)
            cnt = cnt+1;
        end
    end
    if cnt == m
        disp('Unbounded');
        break;
    end
    leave = 1; ratio = 10000000;
    for i=1:m
        if(S(i,enter) >0 && ratio > S(i,2)/S(i,enter))
            leave = i;
            ratio = S(i,2)/S(i,enter);
        end
    end
    piv=S(leave,enter);
    S(leave,2:n+2) = S(leave,2:n+2)/piv;
    S(leave,1) = C(enter-2);
    for j = 1:leave-1
        t = S(j,enter);
        S(j,2:n+2) = S(j,2:n+2)-t*S(leave,2:n+2);
    end
    for j = leave+1:m
        t = S(j,enter);
        S(j,2:n+2) = S(j,2:n+2)-t*S(leave,2:n+2);
    end
end
disp(S)