function [hk] = HT_LMS(par)
%HT_LMS 此处显示有关此函数的摘要
%   此处显示详细说明

%% unpack parameters
lmd=par.lmd;
alpha=par.alpha;
rho=par.rho;
uk=par.uk;
dk=par.dk;

%% iterations start
hk=zeros(size(uk,1),size(uk,2)+1);
hk_aux=zeros(size(uk,1),1);
for kk=1:size(uk,2)
    hk_aux=hk_aux+alpha*uk(:,kk)*(dk(kk)-hk_aux.'*uk(:,kk))-alpha*rho*lmd*(hk_aux-hk(:,kk));
    hk(:,kk+1)=hard_th(hk_aux,sqrt(2/rho));
end
end

function [y]=hard_th(x,th)
y=(abs(x)>=th).*x;
end
