function [nc_correct,sig,data_de] = sdNC_significance_test(nc,struct_res,times,seed)
% Significance test by surrogate data method
    nc_correct=nc;
    state=struct_res.state;
    AR=struct_res.AR;
    Sigma=struct_res.Sigma;
    diminfo_select=struct_res.diminfo_select;
    m=diminfo_select.m;
    n=diminfo_select.n;
    errorflag=struct_res.errorflag;
    sig=zeros(m,m);
    rng(seed);
    ncvaluescell=cell(m,n);
    for i=1:m % effect
        for j=1:m % cause
            if i==j
                continue
            end
            ncvalues=zeros(1,times);
            for k=1:times
                statelen=struct_res.statelen(1,state);
                data_de=zeros(statelen,size(diminfo_select.data,2));
                data_de(1:n,:)=4*(rand(size(diminfo_select.data(1:n,:)))-0.5);
                AR_de=cell(1,n);
                for r=1:n
                    ARone=AR{r,1};
                    ARone(i,j)=0;
                    AR_de{r,1}=ARone;
                end
                for t=n+1:statelen
                    sum1=0;
                    for r=1:n
                        sum1=sum1+(AR_de{r,1}*data_de(t-r,:)')';
                    end
                    data_de(t,:)=sum1+mvnrnd(zeros(1,m),round(Sigma,6),1);
                end
                nc1=NC_data_NHSMM(i,j,data_de(1:size(data_de,1),:),n,errorflag);
                ncvalues(1,k)=nc1;
            end
            ncvaluescell{i,j}=ncvalues;
            pts=abs(nc{j+1,i+1});
            f=ksdensity(ncvalues,pts,'Function','cdf');
            sig(j,i)=double(1-f);
            if (sig(j,i)>=0.05)
                nc_correct{j+1,i+1}=0;
            end
            fprintf("Significance test: finished:%d->%d, state %d\n",j,i,state)
        end
    end        

end

