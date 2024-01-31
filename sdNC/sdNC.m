function [nc,struct_res] = sdNC(struct,diminfo,ABGZE,selectN,state,errorflag)
% calculate casaulity network by state, set errorflag=1 to contain the
% contribution from noise to the effect time course
    n=diminfo.n;
    T=diminfo.T;
    m=diminfo.m;
    data=diminfo.data;
    [~,index]=sort(ABGZE.average_Gamma,'descend');
    indexN=index(1,1:selectN);
    struct_select.W=struct.W(indexN,1);
    struct_select.MN=struct.MN(indexN,1);
    diminfo_select=diminfo;
    diminfo_select.N=selectN;
    gamma=ABGZE.Gamma(:,indexN)./sum(ABGZE.Gamma(:,indexN),2);
    for j=1:n
        AR{j,1}=struct_select.MN{state,1}.M(1+(j-1)*m:m*j,:)';
    end
    Sigma=inv(struct_select.W{state,1}.v*struct_select.W{state,1}.S);
    nc=cell(m+1,m+1);
    nc{1,1}='cause\effect';
    for i=1:m
        nc{1,i+1}=num2str(i);
        nc{i+1,1}=num2str(i);
        nc{i+1,i+1}=0;
    end
    jj=zeros(m,m);
    for i=1:m % effect
        for j=1:m % cause
            sum3=0;
            for t=n+1:T
                sum1=0;
                for r=1:n
                    sum1=sum1+AR{r,1}(i,j)*data(t-r,j);
                end
                sum2=sum1^2;
                sum3=sum3+gamma(t-n,state)*sum2;
            end
            jj(i,j)=sum3;
            fprintf("sdNC finished:%d->%d, state %d\n",j,i,state)
        end
    end
    sum4=sum(jj,2);
    statelen=round(sum(gamma,1));
    if errorflag==1
        for i=1:m
            for j=1:m
                nc{j+1,i+1}=jj(i,j)/(sum4(i)+statelen(1,state)*Sigma(i,i));
            end
        end
    else
        for i=1:m
            for j=1:m
                nc{j+1,i+1}=jj(i,j)/sum4(i);
            end
        end
    end
    struct_res.AR=AR;
    struct_res.Sigma=Sigma;
    struct_res.diminfo_select=diminfo_select;
    struct_res.struct_select=struct_select;
    struct_res.state=state;
    struct_res.errorflag=errorflag;
    struct_res.gamma=gamma;
    struct_res.statelen=statelen;
    struct_res.nctestlen=floor(T/selectN);
end

