function [data,state,data_struct,sample_] = synthetic_data_simuHCP_non(P,T,eachlen,varargin)
% P: transition probability matrix
% T: time points
% eachlen: length of the sample for generating MAR coefficients
% normlevel, dispersionlevel: control the divergence of state-specific MAR coefficients
    vars=["N","D","n_lag",'seed','normlevel','dispersionlevel'];
    values=[5,20,4,2023,2,20];
    p=inputParser;
    for i=1:length(vars)
            addParameter(p,vars(i),values(i));
    end
    parse(p,varargin{:});
    N=p.Results.N;
    D=p.Results.D;
    n_lag=p.Results.n_lag;
    seed=p.Results.seed;
    normlevel=p.Results.normlevel;
    dispersionlevel=p.Results.dispersionlevel;

    load data_group.mat

    dismin=0;
    it=1;
    rng(seed)
    totalN=size(data_group,1)/eachlen;
    while dismin <= dispersionlevel
        sample_=randi(totalN,[N,1]);
        for i=1:N
            index=sample_(i,1);
            data_AR{i,1}=data_group(1+(index-1)*eachlen:index*eachlen,:);

            VAR=gc_my(data_AR{i,1},n_lag);
            AR_cell{i,1}=[];
            AR_cellT{i,1}=[];
            for k=1:n_lag
                AR_cell{i,1}=[AR_cell{i,1};VAR.AR{1,k}];
                AR_cellT{i,1}=[AR_cellT{i,1};VAR.AR{1,k}'];
            end
            var_cell{i,1}=VAR.var;
        end
        normdif_matrix = matrix_norms_dif(AR_cellT,normlevel);
        dismin=min(normdif_matrix,[],'all');
        it=it+1;
        fprintf("it=%d,dismin=%d\n",it,dismin)
    end

   Walk=[];
    dumatrix=zeros(N,D);
    for i=1:N
        A=rand(1,D);
        dumatrix(i,:)=A./sum(A);
    end
    TT=0;
    k=1;
    while TT<=T
        if k==1
            rand_N=randi(N);
            rand_d=rand();
            d=find((cumsum(dumatrix(rand_N,:))-rand_d)>=0);
            d=d(1,1);
            Walk=[Walk;kron(rand_N,ones(d,1))];
            TT=TT+d;
            k=k+1;
        else
            old_N=Walk(end,1);
            runi=rand();
            a=find((cumsum(P(old_N,:))-runi)>=0);
            new_N=a(1,1);
            rand_d=rand();
            d=find((cumsum(dumatrix(new_N,:))-rand_d)>=0);
            d=d(1,1);
            Walk=[Walk;kron(new_N,ones(d,1))];
            TT=TT+d;
        end
    end
    state_seqs=Walk(1:T,1);
    for i=1:N
        state_gamma(i,1)=sum(state_seqs==i)/T;
    end
    data=zeros(T,size(data_group,2));
    for k=n_lag:size(state_seqs,1)
        if k==n_lag
            data(1:n_lag,:)=data_AR{state_seqs(1,1),1}(1:n_lag,:);
        else
            Zt=[];
            for i=1:n_lag
                Zt=[Zt;data(k-i,:)'];
            end
            data(k,:)=(AR_cellT{state_seqs(k,1),1}'*Zt+mvnrnd(zeros(size(data_group,2),1),var_cell{state_seqs(k,1),1},1)')';
        end
    end
    data_struct.data=data_AR;
    data_struct.AR_for_understand=AR_cell;
    data_struct.ART_for_cal=AR_cellT;
    data_struct.var=var_cell;
    state.state_seqs=state_seqs;
    state.state_gamma=state_gamma;
    state.dumatrix=dumatrix;
    state.TT=TT;
end


