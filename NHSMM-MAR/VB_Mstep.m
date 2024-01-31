function struct_update = VB_Mstep(diminfo,struct,abgze_cell)
% struct_update: update model parameters \Theta

    struct_update=struct;
    D=diminfo.D;
    n=diminfo.n;
    N=diminfo.N;
    T=diminfo.T;
    Y=diminfo.Y;
    Z=diminfo.Z;
    gamma_cell=abgze_cell.gamma_cell;
    zeta_cell=abgze_cell.zeta_cell;
    eta_cell=abgze_cell.eta_cell;

    gamma_vec_cell=cell(N,1);
    for i=1:N
        diag1=zeros(T-n,1);
        for k=1:T-n
            sum1=numerical_overflow(sum(gamma_cell{k,1}(i,:)));
            diag1(k,1)=sum1;
        end
        gamma_vec_cell{i,1}=diag1;
    end
    % update the Matrix normal Wishart distribution
    for i=1:N
        % to save memory and avoid large-scale matrix operation
        ZGZ=0;
        ZGY=0;
        YGY=0;
        for k=1:T-n
            ZGZ=ZGZ+Z(k,:)'*gamma_vec_cell{i,1}(k,1)*Z(k,:);
            ZGY=ZGY+Z(k,:)'*gamma_vec_cell{i,1}(k,1)*Y(k,:);
            YGY=YGY+Y(k,:)'*gamma_vec_cell{i,1}(k,1)*Y(k,:);
        end
        struct_update.MN{i,1}.Phi=inv(inv(struct.MN{i,1}.Phi)+ZGZ);
        struct_update.MN{i,1}.M=struct_update.MN{i,1}.Phi*(inv(struct.MN{i,1}.Phi)*struct.MN{i,1}.M+ ...
            ZGY);
        struct_update.W{i,1}.S=inv(inv(struct.W{i,1}.S)+YGY- ...
            struct_update.MN{i,1}.M'*inv(struct_update.MN{i,1}.Phi)*struct_update.MN{i,1}.M+ ...
            struct.MN{i,1}.M'*inv(struct.MN{i,1}.Phi)*struct.MN{i,1}.M);
        sum2=0;
        for k=1:T-n
            sum2=sum2+sum(gamma_cell{k,1}(i,:));
        end
        struct_update.W{i,1}.v=struct.W{i,1}.v+sum2;
    end


    % update \pi
    gamma1_vec=reshape(gamma_cell{1,1}',N*D,1);
    for i=1:N
        for d=1:D
            sum3=0;
            if i==N&&d==D
                sum3=0;
            else
                for j=(i-1)*N+d+1:N*D
                    sum3=sum3+gamma1_vec(j,1);
                end
            end
            struct_update.ppi.alpha(i,d)=struct.ppi.alpha(i,d)+gamma_cell{1,1}(i,d);
            struct_update.ppi.beta(i,d)=struct.ppi.beta(i,d)+sum3;
        end
    end
    % update A
    for i=1:N
        for j=1:N
            if i==j
                continue
            else
                sum4=0;
                for k=1:T-n-1
                    sum4=sum4+zeta_cell{k,1}(i,j);
                end
                sum1=0;
                if j==N
                    sum1=0;
                else
                    for r=j+1:N
                        for k=1:T-n-1
                            sum1=sum1+zeta_cell{k,1}(i,r);
                        end
                    end
                end
                struct_update.a.alpha(i,j)=struct.a.alpha(i,j)+sum4;
                struct_update.a.beta(i,j)=struct.a.beta(i,j)+sum1;
            end
        end
    end

    % update du
    for i=1:N
        for d=1:D
            sum2=0;
            for k=1:T-n-1
                sum2=sum2+eta_cell{k,1}(i,d);
            end
            struct_update.pd.alpha(i,d)=struct.pd.alpha(i,d)+sum2;
            sum5=0;
            if d==D
                sum5=0;
            else
                for dd=d+1:D
                    for k=1:T-n-1
                        sum5=sum5+eta_cell{k,1}(i,dd);
                    end
                end
            end
            struct_update.pd.beta(i,d)=struct.pd.beta(i,d)+sum5;
        end
    end
end

