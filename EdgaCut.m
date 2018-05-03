function  OutputA = EdgaCut(A,f_size)
half_f_size = (f_size-1)/2;
sizeA = size(A);
dim = length(sizeA);
if dim == 1
    OutputA = A(half_f_size+1:end-half_f_size);
else
    if dim == 2
    OutputA = A(half_f_size+1:end-half_f_size , half_f_size+1:end-half_f_size);
    else
        if dim == 3
            OutputA = A(half_f_size+1:end-half_f_size , half_f_size+1:end-half_f_size,:);
        else
            if dim == 4
                OutputA = A(half_f_size+1:end-half_f_size , half_f_size+1:end-half_f_size,:,:);
            else
                if dim == 5
                    OutputA = A(half_f_size+1:end-half_f_size , half_f_size+1:end-half_f_size,:,:,:);
                else
                    disp('dim <= 5')
                end
            end
        end
    end
end