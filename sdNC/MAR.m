function result = MAR(data,p)
% Maximum likelihood estimate of MAR
% p: Autoregressive lag order
    row=size(data,1);
    col=size(data,2);
    ytxtt=0;
    xtxtt=0;
    epsiont_times=0;
    XT=ones(row-p,col*p);
    for i=p+1:row
        yt=data(i,:)';
        xt=[];
        for j=1:p
            xt=[xt data(i-j,:)];
        end
        XT(i-p,:)=xt;
        ytxtt=ytxtt+yt*xt;
        xtxtt=xtxtt+xt'*xt;
    end
    coef=ytxtt/xtxtt;
    for i=p+1:row
        epsiont=data(i,:)'-coef*XT(i-p,:)';
        epsiont_times=epsiont_times+epsiont*epsiont';
    end
    var=epsiont_times/(row-p);
    coef_AR=cell(1,p);
    for i=1:p
        coef_AR{1,i}=coef(:,(i-1)*col+1:i*col);
    end
    result=struct();
    result.var=var;
    result.AR=coef_AR;
    result.XT=XT;
    result.P=p;
    sum1=0;
    for i=p+1:row
        sum1=sum1+(data(i,:)'-coef*XT(i-p,:)')'*inv(var)*(data(i,:)'-coef*XT(i-p,:)');
    end
    loglike=-(row-p)*col/2*log(2*pi)+(row-p)/2*log(det(inv(var)))-0.5*sum1;
    AIC=log(det(var))+2*col^2*p/(row-p);
    BIC=log(det(var))+col^2*p*log(row-p)/(row-p);
    HQ=log(det(var))+2*p*col^2*log(log(row-p))/(row-p);
    result.AIC=AIC;
    result.BIC=BIC;
    result.HQ=HQ;
end

