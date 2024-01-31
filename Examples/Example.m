
% This example script provides some basic demonstration of the framework NSHMM-MAR-sdNC utilization.
% Under Examples/, there are data and scripts used for the analysis conducted in the manuscript.

% the complete code in Simulation1. and Simulation2. may be time-consuming,
% one can quickly get to know the NSHMM-MAR-sdNC utilization form "Basic
% demonstration about NSHMM-MAR-sdNC" and code in Simulation3.

% code about competitive model HMM-MAR is referred to 
% Diego Vidaurre, Andrew J. Quinn, Adam P. Baker, David Dupret, Alvaro Tejero-Cantero and Mark W. Woolrich (2016) Spectrally resolved fast transient brain states in electrophysiological data. NeuroImage. Volume 126, Pages 81–95.
% and toolbox HMM-MAR available at https://github.com/OHBA-analysis/HMM-MAR.


%% Basic demonstration about NSHMM-MAR-sdNC
% 1.synthetic data generate --> 2.Initialization --> 3.NHSMM-MAR modeling and inference --> 4.sdNC network constrcution 

% 1.synthetic data 
% take the synthetic data used in Section3.1. "NHSMM-MAR:robust to infer
% state number" as an example, characterizing by a Markov chian and MAR emisiion model where 
% the coefficients are derived from Preprocessed HCP data (see synthetic_data_simu1.m for details).

addpath(genpath('../NHSMM-MAR-sdNC'))

[data,state,data_struct,sample_] = synthetic_data_simuHCP_geo(transition_matrix_geo(3,0.8,0.9,200),250*3*4,120,N=3,seed=200,dispersionlevel=20);

% 2.Initialization

% [struct,diminfo,ABGZE] = NHSMM_MAR_VB_inferN(data,4,20,10,30);
% To given a suitable state number used in Initialization.m, one may use
% NHSMM_MAR_VB_inferN.m. For simplicity we omit this step to save time. 
[inistruct] = Initialization(data,4,3,30,200); 

% 3. NHSMM-MAR modeling and inference
[struct,diminfo,ABGZE] = NHSMM_MAR_VB(data,inistruct,4,20,10,30);

disp(ABGZE.average_Gamma) % 3 valid states (non-trivial state fractional occupancy) are identified

% matching the states by MAR coefficients
N=3;
dif=zeros(N,N);
for i=1:N
    for j=1:N
        dif(i,j)=sum(abs(data_struct.ART_for_cal{i,1}-struct.MN{j,1}.M),"all");
    end
end
indexmatch = zeros(1,N);
indexstate = [1,2,3]; % by ABGZE.average_Gamma
for i=1:N
    [~,minindex] = min(dif(i,:));
    indexmatch(1,i) = indexstate(minindex);
end
ACC = NHSMM_ACC(data_struct,struct,indexmatch);

% 4.sdNC network constrcution
for i=1:N
    [nc,struct_res] = sdNC(struct,diminfo,ABGZE,N,i,1); % sdNC calculate
    nc_cell{i,1} = nc;
    [nc_correct,sig,data_de] = sdNC_significance_test(nc,struct_res,500,200); % significance test by surrogate data method
    sig_cell{i,1} = sig;
    [index,nc_FDR] = NC_FDR(nc,sig,0.001); % FDR correction
    NC_FIX{i,1} = nc_FDR;
end

%% Simulation1. Section3.1. "NHSMM-MAR:robust to infer state number"

% main=pwd;

for i=3:3
%     dir1=fullfile(main,sprintf("N=%d",i));
%     cd(dir1)
    for j=1:50
        seed=j*100;
%         folderdataname=sprintf("seed%d",seed);
%         mkdir(folderdataname)
%         cd(fullfile(dir1,folderdataname))
        [data,state,data_struct,sample_] = synthetic_data_simuHCP_geo(transition_matrix_geo(3,0.8,0.9,seed),250*3*4,120,N=i,seed=seed,dispersionlevel=20);
%         save("data.mat","data","data_struct","sample_")
%         mkdir("HMM-MAR")
%         mkdir("NHSMM-MAR")
        for k= i:10
            options.K=k;
            options.order=4;
            [hmm, Gamma, Xi, vpath] = hmmmar (data,1000*i,options); % competitive model HMM-MAR
            occupy_rate=sum(Gamma,1)/(1000*i-4);
%             cd(fullfile(dir1,folderdataname,"HMM-MAR"))
%             save(sprintf("resHMM_K=%d",k),"hmm","Gamma")
            [inistruct]=Initialization(data,4,i,30,300);
            [struct,diminfo,ABGZE] = NHSMM_MAR_VB(data,inistruct,4,20,k,30);
%             cd(fullfile(dir1,folderdataname,"NHSMM-MAR"))
%             save(sprintf("resNHSMM_K=%d",k),"struct","diminfo","ABGZE")
            fprintf("i=%d,j=%d,k=%d\n",i,j,k)
        end
%         cd(dir1)
    end
end

% The estimation accuracy (evaluated by VN and ACC) derived from above code
% is summarized in Examples/Statistical analysis and visualization/Simualtion1/.

%% Simulation2. Section3.2. "NHSMM-MAR:free from state duration distribution"

% geometric state duration

% main_geo=fullfile(pwd,"geo");
% mkdir(main_geo)
for i=2:6
%     dir1=fullfile(main_geo,sprintf("N=%d",i));
%     mkdir(dir1)
%     cd(dir1)
    for j=1:50
        seed=j*100;
%         folderdataname=sprintf("seed%d",seed);
%         mkdir(folderdataname)
%         cd(fullfile(dir1,folderdataname))
        [data,state,data_struct,sample_]=synthetic_data_simuHCP_geo(transition_matrix_geo(i,0.93,0.95,seed),i*4*250,120,N=i,seed=seed,dispersionlevel=20);
%         save("data.mat","data","data_struct","sample_")
%         mkdir("HMM-MAR")
%         mkdir("NHSMM-MAR")
        options.K=2*i;
        options.order=4;
        [hmm, Gamma, Xi, vpath] = hmmmar (data,1000*i,options);
        occupy_rate=sum(Gamma,1)/(1000*i-4);
%         cd(fullfile(dir1,folderdataname,"HMM-MAR"))
%         save("resHMM","hmm","Gamma","occupy_rate")
        [inistruct]=Initialization(data,4,i,30,300);
        [struct,diminfo,ABGZE] = NHSMM_MAR_VB(data,inistruct,4,20,i+1,30);
%         cd(fullfile(dir1,folderdataname,"NHSMM-MAR"))
%         save("resNHSMM","struct","diminfo","ABGZE")
        fprintf("i=%d,j=%d\n",i,j)
%         cd(dir1)
    end
end

% nonparametric state duration

% main_non=fullfile(pwd,"non");
% mkdir(main_non)
for i=2:6
%     dir1=fullfile(main_non,sprintf("N=%d",i));
%     mkdir(dir1)
%     cd(dir1)
    for j=1:50
        seed=j*100;
%         folderdataname=sprintf("seed%d",seed);
%         mkdir(folderdataname)
%         cd(fullfile(dir1,folderdataname))
        [data,state,data_struct,sample_] = synthetic_data_simuHCP_non(transition_matrix_non(i,seed),i*4*250,120,D=10,N=i,seed=seed,dispersionlevel=20);
%         save("data.mat","data","data_struct","sample_")
%         mkdir("HMM-MAR")
%         mkdir("NHSMM-MAR")
        options.K=2*i;
        options.order=4;
        [hmm, Gamma, Xi, vpath] = hmmmar (data,1000*i,options);
        occupy_rate=sum(Gamma,1)/(1000*i-4);
%         cd(fullfile(dir1,folderdataname,"HMM-MAR"))
%         save("resHMM","hmm","Gamma","occupy_rate")
        [inistruct]=Initialization(data,4,i,30,300);
        [struct,diminfo,ABGZE] = NHSMM_MAR_VB(data,inistruct,4,20,i+1,30);
%         cd(fullfile(dir1,folderdataname,"NHSMM-MAR"))
%         save("resNHSMM","struct","diminfo","ABGZE")
        fprintf("i=%d,j=%d\n",i,j)
%         cd(dir1)
    end
end

% The estimation accuracy (evaluated by VN and ACC) derived from above code
% is summarized in Examples/Statistical analysis and visualization/Simualtion2/.

%% Simulation3. Section3.3. "NHSMM-MAR-sdNC: application to synthetic data"

N=5;
d=3;

[data,data_struct] = synthetic_data_simu3(transition_matrix_non(5,200),T=N*d*250,sigma2=0.5,D=10,seed=200);

[inistruct]=Initialization(data,3,5,100,300);

[struct,diminfo,ABGZE] = NHSMM_MAR_VB(data,inistruct,3,20,10,30);

dif=zeros(N,N);
for i=1:N
    for j=1:N
        dif(i,j)=sum(abs(data_struct.ART_for_cal{i,1}-struct.MN{j,1}.M),"all");
    end
end

for i=1:N
    [nc,struct_res] = sdNC(struct,diminfo,ABGZE,N,i,1); % sdNC calculate
    nc_cell{i,1} = nc;
    [nc_correct,sig,data_de] = sdNC_significance_test(nc,struct_res,500,200); % significance test by surrogate data method
    sig_cell{i,1} = sig;
    [index,nc_FDR] = NC_FDR(nc,sig,0.001); % FDR correction
    NC_FIX{i,1} = nc_FDR;
end

%% Application on HCP rs-fMRI data

load HCPdata.mat

% optimal autoregressive lag order is 4

% maxlag=15;
% AIC=zeros(maxlag,1);
% BIC=zeros(maxlag,1);
% HQ=zeros(maxlag,1);
% for i=1:maxlag
%     r=MAR(data_group,i);
%     AIC(i,1)=r.AIC;
%     BIC(i,1)=r.BIC;
%     HQ(i,1)=r.HQ;
% end

[inistruct]=Initialization(data_group,4,4,50,300);

[struct,diminfo,ABGZE] = NHSMM_MAR_VB(data_group,inistruct,4,20,10,30);

Gamma=ABGZE.Gamma;
average_Gamma=ABGZE.average_Gamma;

% save("HCPdata_result.mat","average_Gamma","inistruct","diminfo","struct","Gamma")
% Results of the estimation is saved in HCPdata_result.mat

% Herein surrogate data method is avoided to raise efficiency, to induce
% sparsity, the estimated sdNC are thresholded at 20% of their maximal
% value in the manuscript.

% load HCPdata_result.mat
% ABGZE.Gamma=Gamma;
% ABGZE.average_Gamma=average_Gamma;

% VN=5
for i=1:5
    [nc,struct_res] = sdNC(struct,diminfo,ABGZE,N,i,1); % sdNC calculate
    nc_cell{i,1} = nc;
%     [nc_correct,sig,data_de] = sdNC_significance_test(nc,struct_res,500,200); % significance test by surrogate data method
%     sig_cell{i,1} = sig;
%     [index,nc_FDR] = NC_FDR(nc,sig,0.001); % FDR correction
%     NC_FIX{i,1} = nc_FDR;
end

