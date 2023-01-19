function [hk] = ZA_LMS(par)
%ZA_LMS 此处显示有关此函数的摘要
%   此处显示详细说明

%% unpack parameters
uk=par.uk;
dk=par.dk;
mu=par.mu;
gamma=par.gamma;

%% iterations start
hk=zeros(size(uk,1),size(uk,2)+1);
for kk=1:size(uk,2)
    e_k=dk(kk)-uk(:,kk).'*hk(:,kk);
    hk(:,kk+1)=hk(:,kk)-mu*gamma*sign(hk(:,kk))+mu*e_k*uk(:,kk);
end
end