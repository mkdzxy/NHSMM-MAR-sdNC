function geomeans = geo_GD(alpha,beta)
    len=length(alpha);
    for i=1:len
        sum1=0;
        if i>1
            for j=1:i-1
                sum1=sum1+psi(beta(j,1))-psi(alpha(j,1)+beta(j,1));
            end
        else
            sum1=0;
        end
        geomeans(i,1)=exp(psi(alpha(i,1))-psi(alpha(i,1)+beta(i,1))+sum1);
    end
end

