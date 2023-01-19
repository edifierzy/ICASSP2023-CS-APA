function [u] = generate_AR(gamma,K)
%GENERATE_AR 此处显示有关此函数的摘要
%   此处显示详细说明
u=zeros(K,1);
u(1)=randn;
for ii=2:K
    u(ii)=gamma*u(ii-1)+randn;
end
u=u/sqrt(var(u));
end

