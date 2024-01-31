function [abgx_cell,ck] = E_step(diminfo,struct,bik)
    T=diminfo.T;
    N=diminfo.N;
    n=diminfo.n;
%     m=diminfo.m;
%     data=diminfo.data;
    alpha_cell=zeros(N,T-n);
    ck_alpha=zeros(T-n,1);
    beta_cell=zeros(N,T-n);
    ck_beta=zeros(T-n,1);
    gamma_cell=zeros(N,T-n);
    ck_gamma=zeros(T-n,1);
    xi_cell=cell(T-1-n,1);
    ck_xi=zeros(T-1-n,1);
    % alpha
    for k=1:T-n
        if k==1
            for i=1:N
                F(i,1)=numerical_overflow(struct.pi(i,1)*bik(i,1));
            end
            ck_alpha(k,1)=sum(F);
            alpha_cell(:,k)=F/ck_alpha(k,1);
        else
            for i=1:N
                sum1=0;
                for j=1:N
                    sum1=sum1+numerical_overflow(alpha_cell(j,k-1)*struct.a(j,i)*bik(i,k));
                end
                F(i,1)=sum1;
            end
            ck_alpha(k,1)=sum(F);
            alpha_cell(:,k)=F/ck_alpha(k,1);
        end
    end
    % beta
    for k=T-n:-1:1
        if k==T-n
            beta_cell(:,k)=1/N;
            ck_beta(k,1)=sum(beta_cell(:,k));
        else
           for i=1:N
               sum2=0;
               for j=1:N
                    sum2=sum2+numerical_overflow(beta_cell(j,k+1)*struct.a(i,j)*bik(j,k+1));
               end
               F(i,1)=sum2;
           end
           ck_beta(k,1)=sum(F);
           beta_cell(:,k)=F/ck_beta(k,1);
        end
    end
    % gamma
    for k=1:T-n
        for i=1:N
            gamma_cell(i,k)=alpha_cell(i,k)*beta_cell(i,k);
        end
        ck_gamma(k,1)=sum(gamma_cell(:,k));
        gamma_cell(:,k)=gamma_cell(:,k)/ck_gamma(k,1);
    end
    % xi
    for k=1:T-1-n
        for i=1:N
            for j=1:N
                matrix(i,j)=numerical_overflow(beta_cell(j,k+1)*bik(j,k+1)*struct.a(i,j)*alpha_cell(i,k));
            end
        end
        ck_xi(k,1)=sum(matrix,"all");
        for i=1:N
            matrix(i,:)=matrix(i,:)/sum(matrix(i,:));
        end
        xi_cell{k,1}=matrix;
    end
    abgx_cell.alpha=alpha_cell;
    abgx_cell.beta=beta_cell;
    abgx_cell.gamma=gamma_cell;
    abgx_cell.xi=xi_cell;
    ck.ck_alpha=ck_alpha;
    ck.ck_beta=ck_beta;
    ck.ck_gamma=ck_gamma;
    ck.ck_xi=ck_xi;
    abgx_cell.average_gamma=sum(gamma_cell,2)/(T-n);
end

