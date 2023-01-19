function [hk] = l0_NLMS(par)
%L0_NLMS 此处显示有关此函数的摘要
%   此处显示详细说明

%% unpack parameters
uk=par.uk;
dk=par.dk;

gamma=par.gamma; % weight of l0-norm
mu=par.mu; % stepsize
beta=par.beta; % parameter of f_beta

delta=0.01;

%% iterations start
hk=zeros(size(uk,1),size(uk,2)+1);
for kk=1:size(uk,2)
    e_k=dk(kk)-uk(:,kk).'*hk(:,kk);
    hk(:,kk+1)=hk(:,kk)-mu*gamma*f_beta(hk(:,kk),beta)+mu*e_k*uk(:,kk)/(delta+norm(uk(:,kk))^2);
end
end

function [y]=f_beta(x,beta)
y=(abs(x)<=1/beta).*(beta^2*x+sign(x).*beta);
end