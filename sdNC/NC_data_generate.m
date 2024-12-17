function [data,data_struct] = NC_data_generate(varargin)
        vars=["d","trans","seed","sigma2"];
        values=[20,100,1,1];
        p=inputParser;
        for i=1:length(vars)
            addParameter(p,vars(i),values(i));
        end
        parse(p,varargin{:});
        d=p.Results.d;
        trans=p.Results.trans;
        seed=p.Results.seed;
        rng(seed)
        sigma2=p.Results.sigma2;
        transP=randmatrix_N(3,seed);
        Walk=zeros(trans,1);
        for k=1:trans
            if k==1
                r=randi(3);
                Walk(k,1)=r;
            end
            runi=rand(1);
            a=find((cumsum(transP(r,:))-runi)>=0);
            r=a(1,1);
            Walk(k,1)=r;
        end
        state_seqs=kron(Walk,ones(d,1));
        for i=1:3
            state_gamma(i,1)=sum(state_seqs==i)/trans/d;
        end
        AR_cell{1,1}=[0.45,0,0,0,0;
            0.5,0,0,0,0;
            -0.4,0,0,0,0;
            -0.5,0,0,0.35,0.35;
            0,0,0,-0.35,0.35];
        AR_cell{2,1}=[0,0.5,0,0,0;
            0,0.45,0,0,0;
            0,-0.5,0.35,0,0.35;
            0,-0.4,0,0,0;
            0,0,-0.35,0,0.35];
        AR_cell{3,1}=[0,0,0.5,0,0;
            0,0.35,-0.5,0.35,0;
            0,0,0.45,0,0;
            0,-0.35,0,0.35,0;
            0,0,-0.4,0,0];
        for k=1:size(state_seqs,1)
            if k==1
                data(1,:)=mvnrnd(zeros(1,5),diag(sigma2*ones(5,1)));
            else
                data(k,:)=(AR_cell{state_seqs(k,1),1}*data(k-1,:)'+mvnrnd(zeros(5,1),diag(sigma2*ones(5,1)))')';
            end
        end
        data_struct.AR=AR_cell;
        data_struct.transP=transP;
        data_struct.state_seqs=state_seqs;
        data_struct.state_gamma=state_gamma;
        ART_for_cal{1,1}=AR_cell{1,1}';
        ART_for_cal{2,1}=AR_cell{2,1}';
        ART_for_cal{3,1}=AR_cell{3,1}';
        data_struct.ART_for_cal=ART_for_cal;
end

