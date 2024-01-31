function value = numerical_overflow(raw)

    for i=1:size(raw,1)
        for j=1:size(raw,2)
            if raw(i,j)>=realmin('double')
                value(i,j)=raw(i,j);
            else
                value(i,j)=realmin('double');
            end
        end
    end
end

