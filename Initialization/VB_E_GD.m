function [abgze_cell,ck_abgze] = VB_E_GD(diminfo,s)
    D=diminfo.D;
    n=diminfo.n;
    N=diminfo.N;
    T=diminfo.T;
    [alpha_cell,ck_alpha] = alpha_iterate_GD(diminfo,s);
    [beta_cell,ck_beta] = beta_iterate_GD(diminfo,s);
    [gamma_cell,ck_gamma] = gamma_cal_GD(diminfo,alpha_cell,beta_cell);
    [eta_cell,ck_eta] = eta_cal_GD(diminfo,alpha_cell,beta_cell,s);
    [zeta_cell,ck_zeta]= zeta_cal_GD(diminfo,alpha_cell,beta_cell,s);
    abgze_cell.alpha_cell=alpha_cell;
    abgze_cell.beta_cell=beta_cell;
    abgze_cell.gamma_cell=gamma_cell;
    abgze_cell.eta_cell=eta_cell;
    abgze_cell.zeta_cell=zeta_cell;
    ck_abgze.ck_alpha=ck_alpha;
    ck_abgze.ck_beta=ck_beta;
    ck_abgze.ck_gamma=ck_gamma;
    ck_abgze.ck_eta=ck_eta;
    ck_abgze.ck_zeta=ck_zeta;
%     gamma_matrix_cell=cell(N,1);
    % 获取矩阵Gamma_i
%     for i=1:N
%         diag1=zeros(T-n,1);
%         for k=1:T-n
%             sum1=numerical_overflow(sum(gamma_cell{k,1}(i,:)));
%             diag1(k,1)=sum1;
%         end
%         gamma_matrix_cell{i,1}=diag(diag1);
%     end

    Gamma=zeros(T-n,N);
    for i=1:N
        diag1=zeros(T-n,1);
        for k=1:T-n
            sum1=numerical_overflow(sum(gamma_cell{k,1}(i,:)));
            diag1(k,1)=sum1;
        end
        Gamma(:,i)=diag1;
    end
%     for i=1:N
%         Gamma(:,i)=diag(gamma_matrix_cell{i,1});
%     end
    abgze_cell.Gamma=Gamma;
    abgze_cell.average_Gamma=sum(Gamma,1)/(T-n);
end

