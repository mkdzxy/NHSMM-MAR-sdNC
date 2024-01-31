function [struct,diminfo,abgze_cell] = NHSMM_MAR_VB_inferN(data,n,D,N,times,varargin)
% data: T\times m
% inistruct: prior information from initialization
% n: Autoregressive lag order
% D: truncation level \bar{D}
% N: truncation level \bar{N}
% times: maximum iterations
% tol: converged threshold
% delta: change rate of ELBO

% struct: output posterior emission parameters
% diminfo: information about data and settings 
% abgze_cell: marginal statistics used in forward-backward message passing

    vars=["tol","seed","fixtimes"];
    values=[1e-5,2024,times];
    p=inputParser;
    for i=1:length(vars)
        addParameter(p,vars(i),values(i));
    end
    parse(p,varargin{:});
    seed=p.Results.seed;
    tol=p.Results.tol;
    fixtimes=p.Results.fixtimes;
    rng(seed)

    L=[];
    delta=[];
    for i=1:times+1
        if i==1
            [struct,diminfo] = NHSMM_MAR_initialization_inferN(data,n,D,N,seed);
            oldstruct=struct;
            fprintf("Initialization, i=%d\n",i)
        elseif i==2
            s = tildes_GD(oldstruct,diminfo); % Geometric means
            [abgze_cell,ck_abgze] = VB_Estep(diminfo,s); % VB-E: forward-backward recursion and marginal statistic calculation
            struct = VB_Mstep(diminfo,oldstruct,abgze_cell); % VB-M: posterior distribution update
            [l,~]=ELBO(diminfo,ck_abgze.ck_alpha,oldstruct,struct); % ELBO computation
            L=[L;l];
            fprintf("i=%d,ELBO=%d\n",i,L(i-1,1))
            oldstruct=struct;
        else 
            s = tildes_GD(oldstruct,diminfo);
            [abgze_cell,ck_abgze] = VB_E_GD(diminfo,s);
            struct = VB_M_GD(diminfo,oldstruct,abgze_cell);
            [l,~]=L_q_GD(diminfo,ck_abgze.ck_alpha,oldstruct,struct);
            L=[L;l];
            delta_new=(L(i-1,1)-L(i-2,1))/abs(L(i-2,1));
            oldstruct=struct;
            if delta_new>=tol
                delta=[delta;delta_new];
                fprintf("i=%d,ELBO=%d,delta=%d\n",i,L(i-1,1),delta(i-2,1))
            elseif delta_new>=0 && delta_new<tol
                delta=[delta;delta_new];
                fprintf("i=%d,ELBO=%d,delta=%d\n",i,L(i-1,1),delta(i-2,1))
                disp("finished,to tol")
                break
            elseif delta_new<0 && delta(size(delta,1),1)>=0
                delta=[delta;delta_new];
                fprintf("i=%d,ELBO=%d,delta=%d\n",i,L(i-1,1),delta(i-2,1))
                disp("First negative.")
            elseif delta_new<0 && delta(size(delta,1),1)<0
                delta=[delta;delta_new];
                fprintf("i=%d,ELBO=%d,delta=%d\n",i,L(i-1,1),delta(i-2,1))
                disp("Something strange has happened, turn to VB inference with fixed number iterations")
                [struct,diminfo,abgze_cell] = NHSMM_MAR_VB_fixtimes(data,inistruct,n,D,N,fixtimes,seed);
                break
            end
        end
    end
end


