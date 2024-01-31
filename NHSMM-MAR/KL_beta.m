function klbeta = KL_beta(pab,qab)
% pab=[p_alpha,p_beta]'
% qab=[q_alpha,q_beta]'
% KL(p||q)
    klbeta=betaln(qab(1,1),qab(2,1))-betaln(pab(1,1),pab(2,1))-(qab(1,1)-pab(1,1))*psi(pab(1,1))- ... 
        (qab(2,1)-pab(2,1))*psi(pab(2,1))+(qab(1,1)-pab(1,1)+qab(2,1)-pab(2,1))*psi(pab(1,1)+pab(2,1));
end

