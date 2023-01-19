function [hk] = CS_APA_r1_rho1(par)
%CS_APA_R1_RHO1 此处显示有关此函数的摘要
%   此处显示详细说明
% the proposed CS-APA with r=1, rho=1 and constant stepsize

%% unpack parameters
lmd=par.lmd; % the weight of the DC regularizer
eta=par.eta; % the threshold of the debiasing component
uk=par.uk; % the input signal
dk=par.dk; % the output signal
muk=par.muk; % the constant stepsize

%% iterations start
delta=0.01; % the regularization parameter
% the soft shrinkage operator
soft_th=@(x,th)sign(x).*max(0,abs(x)-th);
% the derivative of the debiasing function phi_(eta,0)
phi_prime=@(x,eta)sign(x).*(abs(x)>=eta);

hk=zeros(size(uk,1),size(uk,2)+1);
for kk=1:size(uk,2)
    % compute the APA update vector gk
    ek=dk(kk)-uk(:,kk).'*hk(:,kk);
    gk=uk(:,kk)*ek/(norm(uk(:,kk))^2+delta);
    
    % compute the debiasing component ck
    ck=lmd*phi_prime(hk(:,kk),eta);
    
    % update hk
    hk(:,kk+1)=soft_th(hk(:,kk)+muk*gk+muk*ck,lmd*muk);
end
end

