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

f=mvksdensity(X,xi,'Bandwidth',0.1);




A=1/13*(A+eye(n))-eye(n); %% Mapping eigenvalues in [-1,1]




%% save the gamma coefficients at t_i's
n=900;

t=[-1:2/2^8:1];
nvec=50;
m=length(t);
b=0.01;
M=70; %% Maximum degree of the truncation
r=zeros(m,M+1); %% gamma coefficients at equally spaced points
p=zeros(m,M+1); %% auxillary sequences to find the gamma coefficients


for i=1:m %% t_i's should be given first
    r(i,1)=1/2*(erf((1-t(i))/sqrt(2)*1/b)+erf((1+t(i))/sqrt(2)*1/b));
    r(i,2)=b/sqrt(2*pi)*(-exp(-1/2*(1-t(i))^2*1/b^2)+exp(-1/2*(1+t(i))^2*1/b^2))+t(i)*r(i,1);
    p(i,1)=0;
    p(i,2)=r(i,1);
    for k=2:M %% M degree truncation of Legendre expansion
        if mod(k,2)==0
            r(i,k+1)=1/sqrt(2*pi)*(2*k-1)/k*(-b*exp(-(1-t(i))^2/(2*b^2))-b*exp(-(1+t(i))^2/(2*b^2)))+(2*k-1)/k*b^2*p(i,k)+(2*k-1)/k*t(i)*r(i,k)-(k-1)/k*r(i,k-1);
        else
            r(i,k+1)=1/sqrt(2*pi)*(2*k-1)/k*(-b*exp(-(1-t(i))^2/(2*b^2))+b*exp(-(1+t(i))^2/(2*b^2)))+(2*k-1)/k*b^2*p(i,k)+(2*k-1)/k*t(i)*r(i,k)-(k-1)/k*r(i,k-1);
        end
        p(i,k+1)=(2*k-1)*r(i,k)+p(i,k-1);
    end
end

L=zeros(M+1,nvec);

for i=1:nvec
    v=R(n); %% select random vector
    u=v; 
    w=A*v;
    L(1,i)=v.'*v;  %% save the first two summands for the traces of evaluations of 0 and 1 degree Legendre polynomials
    L(2,i)=v.'*w;
   for j=1:M-1  %% Maximum degree of the Legendre polynomial
       x=(2*j+1)/(j+1)*A*w-j/(j+1)*u;  
       L(j+2,i)=v.'*x; %% save for the trace of evaluation of next degree Legendre polynomial
       u=w; 
       w=x;
   end   
end


l=zeros(M+1,1); %% i-th entry approximates the trace of L_(i-1)(A)

for j=1:M+1
    l(j)=(j-1/2)*sum(L(j,:)); 
end

l=1/nvec*l;


phi=1/n*r*l; %% i-th entry gives the evaluation of phi at t_i, 1/n comes from the definition of DoS

plot(xi,f)
hold on
plot(13*(t+1)-1,phi/13,'--')
hold off

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

