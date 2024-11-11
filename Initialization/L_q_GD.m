function [L,klre] = L_q_GD(diminfo,ck,p,q)
    sum=0;
    for i=1:size(ck,1)
        sum=sum+log(ck(i,1));
    end
    [klb,klbs]=KL_beta_GD(diminfo,p,q);
    mnw=KL_MNW_GD(diminfo,p,q);
    L=-klb-mnw+sum;
    klre.sumlogck=sum;
    klre.mnw=mnw;
    klre.klbs=klbs;
end

