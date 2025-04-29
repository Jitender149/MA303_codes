% why do we need ILPP see thses answers
% write a generalised code for this
A1 = input('Enter matrix for <= type constarint: ');
A2 = input('Enter matrix for = type constraints: ');
A3 = input('Enter matrix for >= type constraint: ');
b = input('Enter the constant term vector for the constraints:');
c = input('Enter coefficients of objective function: ');
x = input('Enter number of varaibles:');
u = x;

[m1,n1] = size(A1); [m2,n2] = size(A2); [m3,n3] = size(A3);
A = zeros(m1+m2+m3,x+m1+m2+2*m3);
if m1 ~= 0
    for i = 1:m1
        A(i,1:n1) = A1(i,:);
    end
end
if m2 ~= 0
    for i = m1+1:m1+m2
        A(i,1:n2) = A2(i-m1,:);
    end
end
if m3 ~= 0
    for i = m1+m2+1:m1+m2+m3
        A(i,1:n3) = A3(i-m1-m2,:);
    end
end
if m1 ~= 0
    for i=1:m1
        A(i,i+n1)=1;
    end
end
if m2 ~= 0
    for i = m1+1:m1+m2
        A(i,i+n2) = 1;
    end
end
if m3 ~= 0
    for i = m1+m2+1:m1+m2+m3
        A(i,i+n3)=1;
    end
end
if m3 ~= 0
    for i = m1+m2+1:m1+m2+m3
        A(i,i+n3+m3) = -1;
    end
end
C = zeros(x+m1+m2+2*m3,1);
for i = 1:x
    C(i) = c(i);
end
if m2 ~= 0
    for i = x+m1+1:x+m1+m2
        C(i) = -5000;
    end
end
if m3 ~= 0
    for i = x+m1+m2+1:x+m1+m2+m3
        C(i) = -5000;
    end
end
disp(A);
disp(C);
[m,n] = size(A);
S = zeros(m,n+2);
for i = m1+1:m1+m2
    S(i,1) = -5000;
end
for i= m1+m2+1:m1+m2+m3
    S(i,1) = -5000;
end
for i = 1:m
    S(i,2) = b(i);
end
for i = 1:m
    for j = 3:n+2
        S(i,j) = A(i,j-2);
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
    x=S(leave,enter);
    S(leave,2:n+2) = S(leave,2:n+2)/x;
    S(leave,1) = C(enter-2);
    for j = 1:leave-1
        t = S(j,enter);
        S(j,2:n+2) = S(j,2:n+2)-t*S(leave,2:n+2);
    end
    for j = leave+1:m
        t = S(j,enter);
        S(j,2:n+2) = S(j,2:n+2)-t*S(leave,2:n+2);
    end
    disp(S);
end
%H = zeros(m,m);
%for i = 1:m1
%   H(:,i) = S(:,i+u+2);
%end
%for i = m1+1:m1+m2
%   H(:,i) = S(:,i+u+2);
%end
%for i = m1+m2+1:m1+m2+m3
%    H(:,i) = S(:,i+u+2);
%end
%disp(H);
%H_ = inv(H);
%H_ = round(H_);
%disp(H_);
%ind= [];
%for i = 1:n
%    for j = 1:m
%       if isequal(A(:,i),H_(:,j))
%          ind = [ind,i];
%         break;
%    end
%end
%end
%disp('The Basic Varaibles Involved:');
%disp(ind);

while true
    [e,d] = size(S);
    maxi = -1; idx = 1;
    for i = 1:e
        if (S(i,2)-floor(S(i,2)) > maxi && i<=2)
            maxi = S(i,2) - floor(S(i,2));
            idx = i;
        end
    end
    if maxi ~= 0
        e = e+1;
        d = d+1;
        S(e,1) = 0;
        S(e,2) = -maxi;
        for i = 3:d-1
            S(e,i) = -(S(idx,i) - floor(S(idx,i)));
        end
        S(e,d) = 1;
        S(1:e-1,d) = 0;
        C(d-2) = 0;
        while true
            disp(S)
            mini = S(1,2); idx = 1;
            for i = 1:e
                if S(i,2) < mini
                    mini = S(i,2);
                    idx = i;
                end
            end
            if mini >=0
                break;
            end
            rat = -11111111111; piv = 3;
            for j = 3:d
                Z_j = 0;
                for i = 1:e
                    Z_j = Z_j + S(i,1)*S(i,j);
                end
                Z_j = Z_j - C(j-2);
                if S(idx,j) < 0 && Z_j/S(idx,j) > rat
                    rat = Z_j/S(idx,j);
                    piv = j;
                end
            end
            S(idx,1) = C(piv-2);
            ele = S(idx,piv);
            S(idx,2:d) = S(idx,2:d)/ele;
            for i = 1:e
                if i ~= idx
                    elem = S(i,piv);
                    S(i,2:d) = S(i,2:d) - elem*S(idx,2:d);
                end
            end
        end
    else
        disp('Optimal Found');
        ansi = 0;
        for i = 1:e
            ansi = ansi + S(i,1)*S(i,2);
        end
        disp(ansi);
        break;
    end
end