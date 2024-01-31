function [alpha_cell,ck_alpha] = alpha_iterate_GD(diminfo,s)
        D=diminfo.D;
        n=diminfo.n;
        N=diminfo.N;
        T=diminfo.T;
        tilde_pi=s.tilde_pi;
        tilde_a=s.tilde_a;
        tilde_b=s.tilde_b;
        tilde_pd=s.tilde_pd;
        alpha_cell=cell(T-n,1);
        ck_alpha=zeros(T-n,1);
        for k=1:T-n
            if k==1
                matrix=zeros(N,D);
                for i=1:N
                    for d=1:D
                      matrix(i,d)=tilde_pi(i,d)*tilde_b(i,1);  
                    end
                end
                ck_alpha(k,1)=sum(matrix,"all");
                matrix=matrix/ck_alpha(k,1);
                alpha_cell{1,1}=matrix;
            else
                matrix=zeros(N,D);
                for i=1:N
                    s=0;
                    for j=setdiff(1:N,i)
                        s=numerical_overflow(s+alpha_cell{k-1,1}(j,1)*tilde_a(j,i));
                    end
                    for d=1:D
                        if(d<D)
                            matrix(i,d)=numerical_overflow(alpha_cell{k-1,1}(i,d+1)*tilde_b(i,k)+s*tilde_b(i,k)*tilde_pd(i,d));
                        end
                        if(d==D)
                            matrix(i,d)=numerical_overflow(s*tilde_b(i,k)*tilde_pd(i,d));
                        end
                    end
                end
                ck_alpha(k,1)=sum(matrix,"all");
                matrix=matrix/ck_alpha(k,1);
                alpha_cell{k,1}=matrix;
            end
        end
end

