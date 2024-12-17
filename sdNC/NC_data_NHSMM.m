function [nc,jj] = NC_data_NHSMM(pointi,pointj,data,p,errorflag)
    result=gc_my(data,p);
    T=size(data,1);
    col=size(data,2);
    nc=cell(col+1,col+1);
    nc{1,1}='cause\effect';
%     nc{col+2,1}='误差';
%     nc_wucha=zeros(1,col);
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
%         for n=1:col
%             if(fuhao(1,n)==0)
%                 nc{n+1,i+1}=0;
%             end
%             if(fuhao(1,n)<0)
%                 nc{n+1,i+1}=-jj(1,n)/(sum(jj)+(T-p)*result.var(i,i));
%             end
%             if(fuhao(1,n)>0)
%                 nc{n+1,i+1}=jj(1,n)/(sum(jj)+(T-p)*result.var(i,i));
%             end
%         end
%     jj(pointi)=0;
    sumjj=sum(jj);
    if errorflag== 1
        nc=jj(pointj)/(sumjj+(T-p)*result.var(pointi,pointi));
    else
        nc=jj(pointj)/sumjj;
    end
%     for i=1:col   
%         nc{col+2,i+1}=nc_wucha(1,i);
%     end
end


