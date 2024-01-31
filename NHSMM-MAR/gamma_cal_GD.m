function [gamma_cell,ck_gamma] = gamma_cal_GD(diminfo,alpha_cell,beta_cell)
    D=diminfo.D;
    n=diminfo.n;
    N=diminfo.N;
    T=diminfo.T;
    gamma_cell=cell(T-n,1);
    ck_gamma=zeros(T-n,1);
    for k=1:T-n
        matrix=zeros(N,D);
        for i=1:N
            for d=1:D
                matrix(i,d)=numerical_overflow(alpha_cell{k,1}(i,d)*beta_cell{k,1}(i,d));
            end
        end
        ck_gamma(k,1)=sum(matrix,"all");
        matrix=matrix/ck_gamma(k,1);
        gamma_cell{k,1}=matrix;
    end
end

