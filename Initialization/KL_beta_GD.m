function [kl_beta,kl] = KL_beta_GD(diminfo,p,q)
%     n=diminfo.n;
    N=diminfo.N;
%     m=diminfo.m;
    D=diminfo.D;
    % KL_pi
    sum1=0;
    for i=1:N
        for j=1:D
            pab=[p.ppi.alpha(i,j),p.ppi.beta(i,j)]';
            qab=[q.ppi.alpha(i,j),q.ppi.beta(i,j)]';
            sum1=sum1+KL_beta(qab,pab);
        end
    end
    KL_pi=sum1;
    % KL_A
    sum2=0;
    for i=1:N
        for j=1:N
            if i==j
                continue
            else
                pab=[p.a.alpha(i,j),p.a.beta(i,j)]';
                qab=[q.a.alpha(i,j),q.a.beta(i,j)]';
                sum2=sum2+KL_beta(qab,pab);
            end
        end
    end
    KL_A=sum2;
    % KL_pd
    sum3=0;
    for i=1:N
        for j=1:D
            pab=[p.pd.alpha(i,j),p.pd.beta(i,j)]';
            qab=[q.pd.alpha(i,j),q.pd.beta(i,j)]';
            sum3=sum3+KL_beta(qab,pab);
        end
    end
    KL_pd=sum3;
    kl_beta=KL_pi+KL_A+KL_pd;
    kl.pi=KL_pi;
    kl.A=KL_A;
    kl.pd=KL_pd;

end

