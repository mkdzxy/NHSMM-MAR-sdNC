function b = bik(struct,diminfo)
    N=diminfo.N;
    T=diminfo.T;
    m=diminfo.m;
    n=diminfo.n;
    Z=diminfo.Z;
    Y=diminfo.Y;
    b=zeros(N,T-n);
    for i=1:N
      for k=1:T-n
          b(i,k)=numerical_overflow((2*pi)^(-m/2)*det(struct.Theta{i, 1}.Sigma)^(-1/2)*exp(-1/2*(Y(k,:)'-struct.Theta{i, 1}.A'*Z(k,:)')'*inv(struct.Theta{i, 1}.Sigma)*(Y(k,:)'-struct.Theta{i, 1}.A'*Z(k,:)')));
      end
    end
    b=b./sum(b,1);
end

