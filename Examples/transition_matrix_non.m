function B = transition_matrix_non(N,seed)
    rng(seed)
    B=zeros(N,N);
    for i=1:N
        A=rand(1,N);
        A(1,i)=0;
        A=A./sum(A);
        for j=1:N
            if j==i
                continue
            else
                B(i,j)=A(1,j);
            end
        end
    end
end

