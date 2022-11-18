N=30;

e   = ones(N,1);
spe = spdiags([e -2*e e], -1:1,N,N);
Iz  = speye(N);
A   = kron(Iz,spe)+kron(spe,Iz);

A = -A;

U = A*0;

for i=1:30
    for j=1:30
        n=(i-1)*30+j;
    
        U(n,n)= exp(- .5*((i-4)^2+(j-5)^2) ) + exp(-.5*((j-15)^2+(i-25)^2));
    end
end

A=A+U*20;

X=eig(A); xi=-1:0.1:25;

g=mvksdensity(X,xi,'Bandwidth',0.001);




A=1/13*(A+eye(n))-eye(n); %% Mapping eigenvalues in [-1,1]


Z=zeros(3,100);


    
t=-1:2/200:1;
m=length(t);
y=zeros(1,m);


M=50;
n=900;
nvec=100;
V=zeros(n,M);
T=zeros(M);
E=zeros(M,nvec);
t1=zeros(M,nvec);


for l=1:nvec
    v=randn(n,1);
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
    for i=3:M
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
    T(M,M-1)=w.'*A*v;
    T(M,M)=w.'*A*w;
    [U,S]=eig(T);
    E(:,l)=eig(T);
    t1(:,l)=U(1,:)';
    
end

s=0.01;

t2=-1:2/200:1;

N1=length(t2);

y=zeros(N1,1);

for i=1:N1
    for l=1:nvec
        for k=1:M
            y(i)=y(i)+t1(k,l)^2*1/sqrt(2*pi)*1/s*exp(-1/(2*s^2)*(t2(i)-E(k,l))^2);
        end
    end
end

y=1/nvec*y;


plot(xi,g)
hold on
plot(13*(t+1)-1,y/13,'--') %% scale funtion y 
hold off

function c=d(x)

c=atan(0.001*x);

end
