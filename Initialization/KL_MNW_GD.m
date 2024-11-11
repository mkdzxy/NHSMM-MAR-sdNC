function klmnw = KL_MNW_GD(diminfo,p,q)
    n=diminfo.n;
    N=diminfo.N;
    m=diminfo.m;
    sum1=zeros(N,1);
    for j=1:N
        for i=1:m
            sum1(j,1)=sum1(j,1)+psi((1+q.W{j,1}.v-i)/2);
        end
    end
    sum2=0;
    for i=1:N
        sum2=sum2+p.W{i,1}.v/2*log(det(p.W{i,1}.S)/det(q.W{i,1}.S))+ ...
            multi_gamma_ln(m,p.W{i,1}.v/2)-multi_gamma_ln(m,q.W{i,1}.v/2)+ ...
            (q.W{i,1}.v-p.W{i,1}.v)/2*sum1(i,1)+0.5*q.W{i,1}.v*(trace(q.W{i,1}.S/p.W{i,1}.S)-m)+ ...
            m/2*log(numerical_overflow(det(p.MN{i,1}.Phi))/numerical_overflow(det(q.MN{i,1}.Phi)))+m/2*trace(q.MN{i,1}.Phi/p.MN{i,1}.Phi)+ ...
            0.5*trace(q.W{i,1}.v*q.W{i,1}.S*(q.MN{i,1}.M-p.MN{i,1}.M)'*inv(p.MN{i,1}.Phi)*(q.MN{i,1}.M-p.MN{i,1}.M))- ...
            0.5*m^2*n;
    end
    klmnw=sum2;
end

