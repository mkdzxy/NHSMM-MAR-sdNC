function lnL = L(diminfo,ck)
    n=diminfo.n;
    T=diminfo.T;
    lnL=0;
    for i=1:T-n
        lnL=lnL+log(ck.ck_alpha(i,1));
    end
end

