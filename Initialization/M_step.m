function struct_new = M_step(diminfo,struct,abgx_cell)
    T=diminfo.T;
    N=diminfo.N;
    n=diminfo.n;
    Z=diminfo.Z;
    Y=diminfo.Y;
%     m=diminfo.m;
    struct_new=struct;
%     for i=1:N
%         Gamma_cell{i,1}=abgx_cell.gamma(i,:);
%     end
    for i=1:N
        Gamma_vec_cell{i,1}=abgx_cell.gamma(i,:);
    end
    for i=1:N
        struct_new.pi(i,1)=abgx_cell.gamma(i,1)/sum(abgx_cell.gamma(:,1));
        for j=1:N
%             sum1=0;
            sum2=0;
            for k=1:T-1-n
%                 sum1=sum1+sum(abgx_cell.xi{k,1}(i,:));
                sum2=sum2+abgx_cell.xi{k,1}(i,j);
            end
            struct_new.a(i,j)=sum2/sum(abgx_cell.gamma(i,:));
        end

        ZGZ=0;
        ZGY=0;
        for k=1:T-n
            ZGZ=ZGZ+Z(k,:)'*Gamma_vec_cell{i,1}(1,k)*Z(k,:);
            ZGY=ZGY+Z(k,:)'*Gamma_vec_cell{i,1}(1,k)*Y(k,:);
        end

%         struct_new.Theta{i,1}.A=inv(Z'*Gamma_cell{i,1}*Z)*(Z'*Gamma_cell{i,1}*Y);
        struct_new.Theta{i,1}.A=inv(ZGZ)*ZGY;
        sum3=0;
        for k=1:T-n
            sum3=sum3+abgx_cell.gamma(i,k)*(Y(k,:)'-struct_new.Theta{i,1}.A'*Z(k,:)')*(Y(k,:)'-struct_new.Theta{i,1}.A'*Z(k,:)')';
        end
        struct_new.Theta{i,1}.Sigma=sum3/sum(abgx_cell.gamma(i,:));
    end
    struct_new.pi=struct_new.pi/sum(struct_new.pi);
    for i=1:N
        struct_new.a(i,:)=struct_new.a(i,:)/sum(struct_new.a(i,:));
    end

end

