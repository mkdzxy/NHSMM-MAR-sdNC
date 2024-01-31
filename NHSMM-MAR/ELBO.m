function [L,klre] = ELBO(diminfo,ck,p,q)
% L: ELBO
% klre: KL divergence
    sum=0;
    for i=1:size(ck,1)
        sum=sum+log(ck(i,1));
    end
    [klb,klbs]=KL_GD(diminfo,p,q); % KL divergence of Generalized Dirichlet distribution
    mnw=KL_MNW(diminfo,p,q); % KL divergence of Matrix normal Wishart distribution
    L=-klb-mnw+sum;
    klre.sumlogck=sum;
    klre.mnw=mnw;
    klre.klbs=klbs;
end

