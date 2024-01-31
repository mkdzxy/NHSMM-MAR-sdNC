function pi = tilde_pi(struct,diminfo)
% Geometric mean of initial distrbution
    N=diminfo.N;
    D=diminfo.D;
    pi_vec=numerical_overflow(geo_GD(reshape(struct.ppi.alpha',N*D,1),reshape(struct.ppi.beta',N*D,1)));
    pi=reshape(pi_vec,D,N)';
%     pi=pi./sum(pi,"all");
end

