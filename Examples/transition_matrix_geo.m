function B = transition_matrix_geo(N,a,b,seed)
% N: state number
% a: lower bound of self-transition probability
% b: upper bound of self-transition probability
    rng(seed)
    B=zeros(N,N);
    for i=1:N
        diagv=a+(b-a)*rand();
        B(i,i)=diagv;
        A=rand(1,N);
        A(1,i)=0;
        A=A./sum(A);
        for j=1:N
            if j==i
                continue
            else
               B(i,j)=(1-diagv)*A(1,j);
            end
        end
    end
end

