function A=Perturb(deltacharge,Gamma,L,alpha)

M=length(alpha);
G1=Gamma*deltacharge;

Gdq=G1(1)*ones(alpha(1),1);

for i=2:M
    Gdq=[Gdq;G1(i)*ones(alpha(i),1)];
end


A=Gdq.*L;

A=L\A;

A=1/2*(A+A');


end