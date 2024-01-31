function normdif_matrix = matrix_norms_dif(matrixcell,normlevel)
    len=size(matrixcell,1);
    normdif_matrix=zeros(len,len);
    for i=1:len
        for j=1:len
            normdif_matrix(i,j)=sum(abs(matrixcell{i,1}-matrixcell{j,1}).^normlevel,"all");
        end
    end
    for i=1:len
        normdif_matrix(i,i)=2022;
    end
end

