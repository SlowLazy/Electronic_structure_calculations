E=eig(H,S);
L=chol(S);


% the exact density
xi=-2:0.01:10;
f=mvksdensity(E,xi,'bandwidth',0.1);
plot(xi,f); fs=f; Es=E;

n= size(H,1); %% dim of the prob
m=30;         %% steps taken for Lanczos

nvec=10;                %% # of random vectors
V=zeros(n,m);           %% V_m
T=zeros(m);             %% T_m

% initial density ...

rho = xi*0;
h   = 0.1; %the bandwidth 

for l=1:nvec
    v=randn(n,1);
    v=L\v;
    sh= v'*H*v/(v'*S*v);
    sh=10^10;
    A = H-sh*S;

    % - compute the basis using Lanczos 
    v=v/norm(v,2);
    V(:,1)=v;

    w=A*v;
    a=w.'*v;
    T(1,1)=a;
    f=w-a*v;
    b=norm(f,2);
    T(1,2)=b;
    w=f/b;
    V(:,2)=w;
    for i=3:m
        a=w.'*A*w;
        T(i-1,i-1)=a;
        b=w.'*A*v;
        T(i-1,i-2)=b;
        f=A*w-a*w-b*v;
        b=norm(f,2);
        T(i-1,i)=b;
        v=w;
        w=f/b;
        V(:,i)=w;
    end
    T(m,m-1)=w.'*A*v;
    T(m,m)=w.'*A*w;
    norm(A*w-T(m,m)*w-T(m,m-1)*v)
    % - projected eigenvalue problem
    H1= V'*H*V; S1=V'*S*V;
    H1=(H1+H1')/2; S1=(S1+S1')/2;
    [Z, D]= eig(H1,S1);
    d= diag(D);
    Z1=S1*Z;
    
    for k=1:m
        rho = rho + 1/(sqrt(2*pi)*h) * exp( -(xi-d(k)).^2/(2*h^2) )*Z1(1,k)*Z1(1,k);
    end
    norm(V*S1-S*V,'fro')
end


rho = rho/nvec; 

z = sum(rho)*0.01;
rho= rho/z;

plot(xi,rho,xi,fs,'--')


