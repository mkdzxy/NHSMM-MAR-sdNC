function [index,nc_FDR] = NC_FDR(nc,nc_sig,pvalue)
% FDR correction
    dim=size(nc_sig,1);
    for i=1:dim
        nc_sig(i,i)=1;
    end
    sig_vec=nc_sig(:);
    [sig_vec_sort,index]=sort(sig_vec);
    sig_vec_sort_pfix=sig_vec_sort(sig_vec_sort<=pvalue);
    for k=1:length(sig_vec_sort_pfix)
        if(k*pvalue/dim/(dim-1)<sig_vec_sort_pfix(k))
            break
        end
    end
    index_vec=index(1:k,1);
    index_dim=reshape(1:dim^2,dim,dim);
    index=zeros(dim,dim);
    for i=1:length(index_vec)
        ii=index_vec(i,1);
        [a,b]=find(index_dim==ii);
        index(a,b)=1;
    end
    for i=1:size(index,1)
        for j=1:size(index,1)
            if index(i,j)==0
                nc{i+1,j+1}=0;
            else
                continue
            end
        end
    end
    nc_FDR=nc;

end

