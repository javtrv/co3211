function T = toeplitz(V)
  
    n = length(V);
    largo = n - ((n/2)-0.5);
    T = 0*ones(largo,largo);
    
    for k=1:1:largo
      for j=0:1:largo-k
        T(j+1,k+j) = V(largo+k-1);
       end
    end
    for k=2:1:largo
      for j=0:1:largo-k
          T(j+k,1+j) = V(largo-k+1);
      end
    end
endfunction