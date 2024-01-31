function [zeta_cell,ck_zeta]= zeta_cal_GD(diminfo,alpha_cell,beta_cell,s)
    D=diminfo.D;
    n=diminfo.n;
    N=diminfo.N;
    T=diminfo.T;
%   tilde_pi=s.tilde_pi;
    tilde_a=s.tilde_a;
    tilde_b=s.tilde_b;
    tilde_pd=s.tilde_pd;
    zeta_cell=cell(T-n-1,1);
    ck_zeta=zeros(T-n-1,1);
    for k=2:T-n
        sumj_d=zeros(N,1);
        for j=1:N
            sum1=0;
            for d=1:D
                sum1=numerical_overflow(tilde_pd(j,d)*beta_cell{k,1}(j,d));
            end
            sumj_d(j,1)=sum1;
        end
        matrix=zeros(N,N);
        for i=1:N
            for j=1:N
                if(i==j)
                    matrix(i,j)=0;
                else
                    matrix(i,j)=numerical_overflow(alpha_cell{k-1,1}(i,1)*tilde_a(i,j)*tilde_b(j,k)*sumj_d(j,1));
                end
            end
        end
        ck_zeta(k-1,1)=sum(matrix,"all");
        matrix=matrix/ck_zeta(k-1,1);
        zeta_cell{k-1,1}=matrix;
    end
end

