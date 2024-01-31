function ACC= NHSMM_ACC(data_struct,struct,index)
% NHSMM_Acc : calculate Acc with NHSMM_Acc
% index: real states with HMM states

    staten=size(data_struct.ART_for_cal,1);
    sum1=0;
    for i=1:staten
        sum1=sum1+sum(abs(data_struct.ART_for_cal{i,1}-struct.MN{index(i),1}.M),"all")/sum(abs(data_struct.ART_for_cal{i,1}),"all");
    end
    ACC=1-sum1/staten;
end
