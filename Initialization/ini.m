function [struct,diminfo] = ini(data,N,n,seed)
    [T,m]=size(data);
    diminfo.data=data;
    diminfo.N=N;
    diminfo.T=T;
    diminfo.m=m;
    diminfo.n=n;
    diminfo.Y=data((n+1):T,:);
    Z=zeros(T-n,m*n);
    for t=1:T-n
        for k=1:n
            Z(t,(k-1)*m+1:k*m)=data(t+n-k,:);
        end
%         Z(t,:)=reshape(data(t:t+n-1,:)',1,m*n);%有问题！
    end
    rng(seed)
    diminfo.Z=Z;
    struct.pi=1/N*ones(N,1);
    struct.a=1/N*ones(N,N);
%     struct.pi=rand(N,1);
%     struct.pi=struct.pi/sum(struct.pi);
%     struct.a=randmatrix_geo(N,0.2,0.8,seed);
    statelen=floor(T/N);
    for i=1:N
%         data_sample=data(randsample(T,round(T/N)),:);
        data_sample=data(((i-1)*statelen+1):(i*statelen),:);
        VAR=gc_my(data_sample,n);
        MAR=[];
        for p=1:n
            MAR=[MAR;VAR.AR{1,p}';];
        end
        struct.Theta{i,1}.A=MAR;
%         struct.Theta{i,1}.A=0.5*rand(m*n,m)-0.25;
%         V=gc_my(data+0.2*(rand(size(data))-0.5),n);
%         struct.Theta{i,1}.Sigma=VAR.var;
        struct.Theta{i,1}.Sigma=cov(data);
%         struct.Theta{i,1}.A=normrnd(0,0.1,m*n,m);
        


    end
end


