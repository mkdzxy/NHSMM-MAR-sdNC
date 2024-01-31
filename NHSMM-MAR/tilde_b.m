function b = tilde_b(struct,diminfo)
% Geometric mean of data likelihood
    N=diminfo.N;
    T=diminfo.T;
    n=diminfo.n;
    m=diminfo.m;
    b=zeros(N,T-n);
    for i=1:N
        for k=n+1:T
            sum_psi=0;
            for ii=1:m
                sum_psi=sum_psi+0.5*psi((struct.W{i,1}.v+1-ii)/2);
            end
            value=exp(-m/2*log(2*pi)+sum_psi+0.5*log(det(struct.W{i,1}.S))- ...
                m/2*trace(struct.MN{i,1}.Phi*diminfo.Z(k-n,:)'*diminfo.Z(k-n,:)) ...
                -0.5*struct.W{i,1}.v*(diminfo.Y(k-n,:)'-struct.MN{i,1}.M'*diminfo.Z(k-n,:)' ...
                )'*struct.W{i,1}.S*(diminfo.Y(k-n,:)'-struct.MN{i,1}.M'*diminfo.Z(k-n,:)'));
            b(i,k-n)=numerical_overflow(value);
        end
    end
    b=b./sum(b,1);
end

