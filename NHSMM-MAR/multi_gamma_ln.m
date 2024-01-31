function value = multi_gamma_ln(m,z)
    sum=0;
    for i=1:m
        sum=sum+gammaln(z+(1-i)/2);
    end
    value=m*(m-1)/4*log(pi)+sum;
end

