function [nc,jj] = NC_surrogatedata(pointi,pointj,data,p,errorflag)
% Casaulity network of surrogate data estimated by NC
    result=gc_my(data,p);
    T=size(data,1);
    col=size(data,2);
    nc=cell(col+1,col+1);
    nc{1,1}='effect\cause';
    for i=1:col
        nc{1,i+1}=num2str(i);
        nc{i+1,1}=num2str(i);
        nc{i+1,i+1}=0;
    end
    jj=zeros(1,col);
    for j=1:col
        sumt=0;
        for t=p+1:T
            sumj1=0;
            sumj2=0;
            for m=1:p
                sumj1=sumj1+result.AR{1,m}(pointi,j)*data(t-m,j);
            end
            sumj2=sumj1^2;
            sumt=sumt+sumj2;
        end
        jj(j)=sumt;
    end
    sumjj=sum(jj);
    if errorflag== 1
        nc=jj(pointj)/(sumjj+(T-p)*result.var(pointi,pointi));
    else
        nc=jj(pointj)/sumjj;
    end
end


