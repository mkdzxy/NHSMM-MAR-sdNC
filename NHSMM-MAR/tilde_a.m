function a = tilde_a(struct,diminfo)
% Geometric mean of transition probability matrix
    N=diminfo.N;
    a=zeros(N,N);
    for i=1:N
        jindex=setdiff(1:N,i);
        a(i,jindex)=geo_GD(struct.a.alpha(i,jindex)',struct.a.beta(i,jindex)')';
    end
%     a=a./sum(a,2);
end

