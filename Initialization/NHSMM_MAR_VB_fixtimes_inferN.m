function [struct,diminfo,abgze_cell] = NHSMM_MAR_VB_fixtimes_inferN(data,n,D,N,times,seed)

    for i=1:times+1
        if i==1
            [struct,diminfo] = NHSMM_MAR_initialization_inferN(data,n,D,N,seed);
            oldstruct=struct;
            fprintf("Fix number VB iterations, i=%d\n",i)
        else 
            s = tildes_GD(oldstruct,diminfo);
            [abgze_cell,~] = VB_Estep(diminfo,s);
            struct = VB_Mstep(diminfo,oldstruct,abgze_cell);
            oldstruct=struct;
            fprintf("Fix number VB iterations, i=%d\n",i)
        end
    end
end


