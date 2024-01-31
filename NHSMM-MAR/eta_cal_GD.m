function [eta_cell,ck_eta] = eta_cal_GD(diminfo,alpha_cell,beta_cell,s)
    D=diminfo.D;
    n=diminfo.n;
    N=diminfo.N;
    T=diminfo.T;
%   tilde_pi=s.tilde_pi;
    tilde_a=s.tilde_a;
    tilde_b=s.tilde_b;
    tilde_pd=s.tilde_pd;
    eta_cell=cell(T-n-1,1);
    ck_eta=zeros(T-n-1,1);
    for k=2:T-n
        sumj_i=zeros(N,1);
        for j=1:N
            sum1=0;
            for i=1:N
                sum1=numerical_overflow(sum1+alpha_cell{k-1,1}(i,1)*tilde_a(i,j));
            end
            sumj_i(j,1)=sum1;
        end
        matrix=zeros(N,D);
        for j=1:N
            for d=1:D
                matrix(j,d)=numerical_overflow(sumj_i(j,1)*tilde_b(j,k)*tilde_pd(j,d)*beta_cell{k,1}(j,d));
            end
        end
        ck_eta(k-1,1)=sum(matrix,"all");
        matrix=matrix/ck_eta(k-1,1);
        eta_cell{k-1,1}=matrix;
    end
end

