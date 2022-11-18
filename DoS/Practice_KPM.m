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

f=mvksdensity(X,xi,'Bandwidth',0.05);




A=1/13*(A+eye(n))-eye(n); %% Mapping eigenvalues in [-1,1]


nvec=50; %% number of random vectors
n=900; %% dimension of matrix
M=100; %% degree of the approximation in Chebyshev polynomial
K=zeros(M+1,nvec); 

for i=1:nvec
    v=R(n); %% select random vector
    u=v; 
    w=A*v;
    K(1,i)=v.'*v;  %% save the first two summands for the traces of 0 and 1 degrees Chebyshev polynomials
    K(2,i)=v.'*w;
   for j=1:M-1    
       z=2*A*w-u;  
       K(j+2,i)=v.'*z; %% save for the trace of next degree Chebyshev polynomial
       u=w; 
       w=z;
   end   
end


a=zeros(M+1,1); %% the averages approximate the traces of T_k(A)'s
a(1)=1./(n*pi).*1/nvec.*sum(K(1,:));
for i=2:M+1
    a(i)=2./(n*pi).*1/nvec*sum(K(i,:));
end



p=M; %% Jackson damping coefficient
ap= pi/(p+2);
for i=1:p
    g(i)= (1-i/(p+2))*sin(ap)*cos(i*ap)+1/(p+2)*cos(ap)*sin(i*ap);
end
g = g/sin(ap); 



x=1/13*(xi+1)-1; %% linear correspondence on the domain
y=a(1).*1./sqrt(1-x.^2);
y=y+a(2).*x.*g(1).*1./sqrt(1-x.^2);
for i=2:M
    y=y+a(i+1).*s(i,x).*g(i).*1./sqrt(1-x.^2);
end

plot(xi,f)
hold on
plot(13*(x+1)-1,y/13) %% scale funtion y 
hold off



function C=s(k,x)
T0=1;
T1=x;
for i=3:k+1
    T=2.*x.*T1-T0;
    T0=T1;
    T1=T;
end
C=T;
% value of k degree Chebyshev polynomial at x
end





function x=R(n)  %% Rademacher random vector

x=rand(n,1);
for i=1:n
    if x(i)<0.5
        x(i)=-1;
    else
        x(i)=1;
    end
end

end