function bitstring = symbolsToBits(symbols)
    N = size(symbols,2);
    bitstring = zeros(1,N*4);
    for i=1:N
       switch (symbols{i})
           case 'A'
               bitstring((i-1)*4+1:i*4)=[0 0 0 0];
           case 'B'
               bitstring((i-1)*4+1:i*4)=[0 0 0 1];
           case 'C'
               bitstring((i-1)*4+1:i*4)=[0 0 1 0];
           case 'D'
               bitstring((i-1)*4+1:i*4)=[0 0 1 1];
           case 'E'
               bitstring((i-1)*4+1:i*4)=[0 1 0 0];
           case 'F'
               bitstring((i-1)*4+1:i*4)=[0 1 0 1];
           case 'G'
               bitstring((i-1)*4+1:i*4)=[0 1 1 0];
           case 'H'
               bitstring((i-1)*4+1:i*4)=[0 1 1 1];
           case 'I'
               bitstring((i-1)*4+1:i*4)=[1 0 0 0];
           case 'J'
               bitstring((i-1)*4+1:i*4)=[1 0 0 1];
           case 'K'
               bitstring((i-1)*4+1:i*4)=[1 0 1 0];
           case 'L'
               bitstring((i-1)*4+1:i*4)=[1 0 1 1];
           case 'M'
               bitstring((i-1)*4+1:i*4)=[1 1 0 0];
           case 'N'
               bitstring((i-1)*4+1:i*4)=[1 1 0 1];
           case 'O'
               bitstring((i-1)*4+1:i*4)=[1 1 1 0];
           case 'P'
               bitstring((i-1)*4+1:i*4)=[1 1 1 1];
       end
    end

end