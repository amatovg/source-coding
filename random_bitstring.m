    function bitstring = random_bitstring(length)
        % Function that generates a random bitsequence.
        % input:
        % length : the length of the bitsequence
        % output:
        % bitsing: a random generated bitstring
        bitstring = round(rand(1,length));
    end