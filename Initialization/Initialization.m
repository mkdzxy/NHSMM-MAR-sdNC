function [struct,diminfo,lnL,abgx_cell,ck] = Initialization(data,n,N,maxit,seed)
    % Initialize with a AR-HMM learned by EM algorithm
    lnL=zeros(maxit,1);
    for i=1:maxit
        if i==1
            [struct,diminfo] = ini(data,N,n,seed);
            fprintf("Initialized\n")
        else
            b = bik(struct,diminfo);
            [abgx_cell,ck] = E_step(diminfo,struct,b);
            struct = M_step(diminfo,struct,abgx_cell);
            lnL(i,1)=L(diminfo,ck);
            fprintf("i=%d,lnL=%d\n",i,lnL(i,1))
        end
    end
end

