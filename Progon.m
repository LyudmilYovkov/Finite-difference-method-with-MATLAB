function res = Progon(A, d)
% b vector
b = diag(A);
s = size(A);
% a vector
for j = 1 : s(1, 2)-1
    a(j) = A(j+1, j);
end
a = [0, a];
% c vector
for i = 1 : s(1, 1)-1
    c(i) = A(i, i+1);
end
c = [c, 0];
alpha(1) = -c(1)/b(1);
beta(1) = d(1)/b(1);
n = length(d);
for k = 2 : n
    alpha(k) = -c(k)/(b(k)+a(k)*alpha(k-1));
    beta(k) = (d(k)-a(k)*beta(k-1))/(b(k)+a(k)*alpha(k-1));
end
x(n) = beta(n);
for k = n-1 : -1 : 1
    x(k) = alpha(k)*x(k+1)+beta(k);
end
res = x';
end

