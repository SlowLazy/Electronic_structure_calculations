function [n1,V,T]=Lanczos(A,v,m)
T=zeros(m);
n=size(A,1);
V=zeros(n,m);
n1=norm(v);
v=v/n1;
V(:,1)=v;
w=A*v;
a=w'*v;
T(1,1)=a;
f=w-a*v;
b=norm(f);
T(1,2)=b;
w=f/b;
V(:,2)=w;
for j=3:m
    a=w'*A*w;
    T(j-1,j-1)=a;
    f=A*w-a*w-b*v;
    b=norm(f);
    T(j-1,j)=b;
    v=w;
    w=f/b;
    V(:,j)=w;
end
T(m,m)=w'*A*w;
T=T+T'-diag(diag(T));
    
end