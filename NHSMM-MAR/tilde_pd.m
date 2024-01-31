function pd = tilde_pd(struct,diminfo)
% Geometric mean of state duration
    N=diminfo.N;
    D=diminfo.D;
    pd=zeros(N,D);
    for i=1:N
        pd(i,:)=geo_GD(struct.pd.alpha(i,:)',struct.pd.beta(i,:)')';
    end
%     pd=pd./sum(pd,2);
end

