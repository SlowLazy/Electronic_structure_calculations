%
% scc calcution using full diagonalization
clear all

load gamma
load A0
load L


beta=1052.58171981327;
mu=-0.164847981649556;
degree=20;
nvec=1;
Gamma=gammaMat;
d = diag(Gamma);
Gamma=(Gamma+Gamma')-diag(d);
m=2;
iterate=10;

M=size(Gamma,1);
alpha=zeros(M,1);

for i=1:M
    alpha(i)=4;
end
n=4*M;
qbar=4*ones(M,1);
Q=zeros(M,iterate);
Q1=zeros(M,iterate);
X=zeros(M,m);
G=zeros(M,m);


lambda=1;

qIn= qbar;

Q1(:,1)=qIn;
X(:,1)=qIn;

qOut=StoLan(A0,L,beta,mu,nvec,degree,alpha);

Q(:,1)=qOut;
G(:,1)=qOut;

for j=2:m
    
    qmix=(1-1/(100))*qIn+1/(100)*qOut;
    
    deltacharge=qmix-qbar;
    
    A1=Perturb(deltacharge,Gamma,L,alpha);
    A=A0+A1;
    
    qIn = qmix;
    
    Q1(:,j)=qIn;
    X(:,j)=qIn;
    qOut=StoLanuni(A,L,beta,mu,nvec,degree,alpha);
    
    Q(:,j)=qOut;
    G(:,j)=qOut;
    lambda=lambda+1;
end


for j=m+1:iterate
    R=G-X;
    a=pinv(R(:,end)-R(:,1:end-1))*R(:,end);
    a=[a;1-sum(a)];
    
    norm(R*a)
    qmix=(1-1/100)*(X*a)+1/100*(G*a);
    Q1(:,j)=qmix;
    X=[X(:,2:end),qmix];
    
    deltacharge=qmix-qbar;
    A1=Perturb(deltacharge,Gamma,L,alpha);
    A=A0+A1;
    
    qOut=StoLanuni(A,L,beta,mu,nvec,degree,alpha);
    Q(:,j)=qOut;
    G=[G(:,2:end),qOut];
    lambda=lambda+1;
end

load('Q1exact.mat')

a=zeros(iterate,1);

for i=1:iterate
a(i)=norm(Q1exact(:,end)-Q1(:,i),'inf');
end

plot(1:iterate,a)