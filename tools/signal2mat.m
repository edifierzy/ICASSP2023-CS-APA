function [X] = signal2mat(x,N)
%SIGNAL2MAT 此处显示有关此函数的摘要
%   此处显示详细说明
X=zeros(N,length(x)-N+1);
for ii=1:(length(x)-N+1)
    X(:,ii)=x(ii:ii+N-1);
end
end

