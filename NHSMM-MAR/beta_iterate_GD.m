function [beta_cell,ck_beta] = beta_iterate_GD(diminfo,s)
    D=diminfo.D;
    n=diminfo.n;
    N=diminfo.N;
    T=diminfo.T;
%   tilde_pi=s.tilde_pi;
    tilde_a=s.tilde_a;
    tilde_b=s.tilde_b;
    tilde_pd=s.tilde_pd;
    beta_cell=cell(T-n,1);
    ck_beta=zeros(T-n-1,1);
    for k=T-n:-1:1
        if k==T-n
            matrix=1/N/D*ones(N,D);
            beta_cell{k,1}=matrix;
            ck_beta(k,1)=sum(matrix,"all");
        else
            matrix=zeros(N,D);
            sumd_j=zeros(N,1);
            for j=1:N
                sum1=0;
                for d=1:D
                    sum1=numerical_overflow(sum1+tilde_pd(j,d)*beta_cell{k+1,1}(j,d));
                end
                sumd_j(j,1)=sum1;
            end
            for i=1:N
                for d=1:D
                    if(d>1)
                        matrix(i,d)=numerical_overflow(tilde_b(i,k+1)*beta_cell{k+1,1}(i,d-1));
                    end
                    if(d==1)
                        sum2=0;
                        for j=1:N
                            sum2=numerical_overflow(sum2+tilde_a(i,j)*tilde_b(j,k+1)*sumd_j(j,1));
                        end
                        matrix(i,d)=sum2;
                    end
                end
            end
            ck_beta(k,1)=sum(matrix,"all");
            matrix=matrix/ck_beta(k,1);
            beta_cell{k,1}=matrix;
        end
    end
end

