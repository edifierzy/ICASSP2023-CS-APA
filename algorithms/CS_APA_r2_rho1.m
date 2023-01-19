function [hk] = CS_APA_r2_rho1(par)
%CS_APA_R2_RHO0 此处显示有关此函数的摘要
%   此处显示详细说明
% the proposed CS-APA with r=2, rho=1 and variable stepsize

%% unpack parameters
lmd=par.lmd; % the weight of the DC regularizer
eta=par.eta; % the threshold of the debiasing component
% hyperparameters for variable stepsize
C=par.C;
mu_max=par.mu_max;
alpha=par.alpha;

uk=par.uk;
dk=par.dk;

%% iterations start
delta=0.01;
% the soft shrinkage operator
soft_th=@(x,th)sign(x).*max(0,abs(x)-th);
% the derivative of the debiasing function phi_(eta,0)
phi_prime=@(x,eta)sign(x).*(abs(x)>=eta);

hk=zeros(size(uk,1),size(uk,2)+1);
Uk=zeros(size(uk,1),2);
dk_vec=zeros(2,1);
p=zeros(2,1);

for kk=1:size(uk,2)
    % update the Uk and dk_vec
    Uk=[uk(:,kk),Uk(:,1)];    
    dk_vec=[dk(kk);dk_vec(1:end-1)];    
    
    % update the variable stepsize mu
    e_k=dk_vec-Uk.'*hk(:,kk);
    Uk_pinvUkUk_ek=Uk*pinv(Uk.'*Uk+delta*eye(2))*e_k;    
    if kk==1
        p=Uk_pinvUkUk_ek;
    else
        p=alpha*p+(1-alpha)*Uk_pinvUkUk_ek;
    end
    mu=mu_max*norm(p)^2/(C+norm(p)^2);
    
    % update hk
    hk(:,kk+1)=soft_th(hk(:,kk)+mu*(Uk_pinvUkUk_ek+lmd*phi_prime(hk(:,kk),eta)),mu*lmd);    
end
end
