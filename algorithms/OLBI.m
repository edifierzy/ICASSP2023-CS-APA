function [hk] = OLBI(par)
%OLBI 此处显示有关此函数的摘要
%   此处显示详细说明
%% unpack parameters
uk=par.uk;
dk=par.dk;
delta=par.delta;
gamma=par.gamma;

%% iterations start
hk=zeros(size(uk,1),size(uk,2)+1);
mk=zeros(size(uk,1),1);
for kk=1:size(uk,2)
    mk=mk+delta*(dk(kk)-hk(:,kk).'*uk(:,kk))*uk(:,kk);
    hk(:,kk+1)=soft_th(mk,gamma);
end
end

function [output]=soft_th(input, threshold)
output=sign(input).*max(0,abs(input)-threshold);
end
