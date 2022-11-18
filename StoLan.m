function q=StoLan(A,L,beta,mu,nvec,degree,alpha)
M=length(alpha);
n=size(A,1);

V=zeros(n,degree);
V1=zeros(n,nvec); %% save r.v.'s

for i=1:nvec
    V1(:,i)=Ra(n);
end

V2=L\V1; %% save initial vectors

for i=1:nvec
    v=V2(:,i);
    [n1,V,T]=Lanczos(A,v,degree);
    [U,D]=eig(T);   
    d=diag(D); 
    f=2./(1+exp(beta*(d-mu)));
    fhT=U*diag(f)*U';
    V2(:,i)=n1*V*fhT(:,1); %% f(A)inv(L)v_i
 
end

V2=L*V2; %% Lf(A)inv(L)v_i

v=v*0;

for i=1:nvec
    v=v+V1(:,i).*V2(:,i); %% use of the diagonal estimator
end


v=1/nvec*v;
q=zeros(M,1);
q(1)=sum(v(1:alpha(1)));
Index=alpha(1);
for i=2:M
    q(i)=sum(v(Index+1:Index+alpha(i)));
    Index=Index+alpha(i);
end

end


function y=Ra(n)  %% Rademacher random vector

y=rand(n,1);
for i=1:n
    if y(i)<0.5
        y(i)=-1;
    else
        y(i)=1;
    end
end

end