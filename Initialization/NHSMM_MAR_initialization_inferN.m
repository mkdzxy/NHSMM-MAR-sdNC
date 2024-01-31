function [struct,diminfo] = NHSMM_MAR_initialization_inferN(data,n,D,N,seed,varargin)
% Initialized by kmeans and MAR to infer the state number 
% data: T\times m
% inistruct: prior information from initialization
% n: Autoregressive lag order
% D: truncation level \bar{D}
% N: truncation level \bar{N}
% X_alpha, X_\beta: parameters of generalized Dirichlet (GD) distribution,
% X \in {\pi, A, du}
    vars=["pi_beta","A_beta","pd_beta"];
    values=[1,1,1];
    p=inputParser;
    for i=1:length(vars)
        addParameter(p,vars(i),values(i));
    end
    parse(p,varargin{:});
    pi_beta=p.Results.pi_beta;
    A_beta=p.Results.A_beta;
    pd_beta=p.Results.pd_beta;
    rng(seed)
    [T,m]=size(data);
    struct.W=cell(N,1);
    mn=m*n;
    struct.MN=cell(N,1);
    iniindex=kmeans(data,N);
    for i=1:N
        struct.MN{i,1}.M=zeros(mn,m);
        datai=data(iniindex==i,:);
        VAR=gc_my(datai,n);
        struct.W{i,1}.S=inv(cov(datai))/(m+2);
        struct.W{i,1}.v=m+2;
        MAR=[];
        for p=1:n
            MAR=[MAR;VAR.AR{1,p}';];
        end
        struct.MN{i,1}.Phi=eye(mn);
        struct.MN{i,1}.M=MAR;
    end
    struct.ppi.alpha=ones(N,D);
    struct.ppi.beta=pi_beta*ones(N,D);
    struct.a.alpha=ones(N,N);
    struct.a.beta=A_beta*ones(N,N);
    for i=1:N
        struct.a.alpha(i,i)=0;
        struct.a.beta(i,i)=0;
    end
    struct.pd.alpha=ones(N,D);
    struct.pd.beta=pd_beta*ones(N,D);
    diminfo.n=n;
    diminfo.D=D;
    diminfo.N=N;
    diminfo.T=T;
    diminfo.m=m;
    diminfo.data=data;
    diminfo.Y=data((n+1):T,:);
    Z=zeros(T-n,m*n);
    for t=1:T-n
        for k=1:n
            Z(t,(k-1)*m+1:k*m)=data(t+n-k,:);
        end
    end
    diminfo.Z=Z;
end



