classdef Channel_Coding
          
   methods(Static=true)
       
        % Zelf geschreven hulpfuncties
       
       function bitstring = random_bitstring(length)
            % Function that generates a random bitsequence.
            % input:
            % length : the length of the bitsequence
            % output:
            % bitsing: a random generated bitstring
            bitstring = round(rand(1,length));
       end

       function generatorMatrix = Generator_Matrix(generator)
            % Function that generates the systematic form of the generator
            % matrix for the (15,11)Hamming code.
            % input:
            % generator: the generatorpolynomial in its bitform
            % n: the length of the codeword
            % k: the length of the informationword
            % output:
            % generatorMatrix: generatormatrix

            n = 15;
            k = 11;
            % generating the generatormatrix 
            generatorMatrix = zeros(k,n);
            for i=1:k
                generatorMatrix(i,i:i+length(generator)-1) = generator;
            end
       end

       function generatorMatrixSys = Systematic_Generator_Matrix(generator)
            % Function that generates the systematic form of the generator
            % matrix for the (15,11)Hamming code.
            % input:
            % generator: the generatorpolynomial in its bitform
            % output:
            % generatorMatrixSys: generatormatrix in its systematic form

            n = 15;
            k = 11;

            % generating the generatormatrix 

            generatorMatrix = Channel_Coding.Generator_Matrix(generator);

            % transformations of the generatormatrix in its systematical
            % form

            n=15;
            k=11;
            generatorMatrixSys = mod(rref(generatorMatrix),2);

        end

       function checkMatrix = Check_Matrix(check)
            % Function that generates the systematic check matrix for the
            % (15,11) Hamming code.
            % input:
            % check: the checkpolynomial in its bitform (IN REVERSE ORDER!!!)
            % output:
            % checkMatrix: checkmatrix

            n = 15;
            k = 11;

            % generating the systematic checkmatrix
            checkMatrix = zeros(n-k,n);
            for i=1:n-k
                checkMatrix(i,i:i+length(check)-1) = check;
            end
        end

       function checkMatrixSys = Systematic_Check_Matrix(generator)
            % Function that generates the systematic check matrix for the
            % (15,11) Hamming code.
            % input:
            % generator: the generatorpolynomial in its bitform
            % output:
            % checkMatrix: checkmatrix in its systematic form

            n = 15;
            k = 11;

            % generating the systematic generatormatrix
            generatorMatrixSys = Channel_Coding.Systematic_Generator_Matrix(generator);

            % generating the checkmatrix from the systematical
            % generatormatrix
            P = generatorMatrixSys(:,k+1:n);
            checkMatrixSys = [P' eye(n-k)];

       end
       
       function [syndromes,cosets] = Syndrome_Table(Hsys)
            % Function that computes the syndrome table, for our specific case
            % input:
            % Hsys: systematic check matrix of the used code
            % output:
            % syndromes: left part of the table, containing the syndromes
            % cosets: right part of the table, containg the coset-leaders

            v = size(Hsys);
            n = v(2);
            %k = v(2)-v(1);
            lengte = 2^v(1); % v(1) = n-k
            syndromes = zeros(lengte,v(1));
            cosets = zeros(lengte,n);

            for i = 2:(lengte)
                for j = 1:v(1)
                    syndromes(i,j)=Hsys(j,n-i+2);
                    cosets(i,n-i+2)=1;
                end
            end

        % deze functie werkt enkel en alleen doordat n+1=2^k=16 en de minimale
        % hammingafstand 3 is waardoor er geen identieke kolommen in onze
        % checkmatrix zitten.
        end
           
        % Functies voor Hamming codering
        
        function bitenc = Ham_encode(bitstring)
            % Functie die een bitstring encodeerd met de (11,15) Hamming
            % code en generatorveelterm 1+x^3+x^4.
            % input:
            % bitstring: vector van ongecodeerde bits met lengte 1xN
            n = 15;
            k = 11;
            
            bitstring = bitstring(:)';  % zorg ervoor dat de input een rijvector is
            N = length(bitstring);
            N_codewords = ceil(N/k);            
            
            % voeg nullen aan de bitstring toe indien N niet deelbaar is
            % door een geheel aantal codewoorden
            bitstring = [bitstring zeros(1, N_codewords*k-N)];                      
           
            % generatorpolynomial in its bit-form
            generator = [1 0 0 1 1];
            
            % generating the systematic generatormatrix
            generatorMatrixSys = Channel_Coding.Systematic_Generator_Matrix(generator);
                       
            % generating the checkmatrix
            %checkMatrix = Channel_Coding.Systematic_Check_Matrix(generator);              
            
            % output: de geencodeerde bits: lengte 15*N_codewords
            bitenc = zeros(1, n*N_codewords);
            
            % encoding the informationwords of bitstring into bitenco 
            for i=1:N_codewords
                bitenc((i-1)*n+1:i*n) = mod(bitstring((i-1)*k+1:i*k) * generatorMatrixSys,2);
            end                      
            
        end

        function bitdec = Ham_decode(bitenc)
            % Functie die een Hammingcode decodeert (foutcorrectie)
            % input:
            % bitenc: vector met gecodeerde bits: lengte moet deelbaar zijn door 15
           
            
            generator = [1 0 0 1 1]; % bits afkomstig van de generatorveelterm
            % generatorMatrixSys = Channel_Coding.Systematic_Generator_Matrix(generator);
            checkMatrix = Channel_Coding.Systematic_Check_Matrix(generator); 
            [syndromes,cosets]=Channel_Coding.Syndrome_Table(checkMatrix);
            
            bitenc = bitenc(:)';    % zorg ervoor dat de input een rijvector is
            N = length(bitenc);
            N_codewords = N/15;
            
            % output: de gedecodeerde bits: lengte 11*N_codewords
            bitdec = zeros(1, 11*N_codewords); 
            
            if(mod(N, 15) ~= 0)
                error('input is geen geheel aantal codewoorden.');
            end
                        
            n = 15;
            k = 11;
            
            for i=1:N_codewords
                r = bitenc((i-1)*n+1:i*n); % deel de geëncodeerde bitsequentie op in meerdere codewoorden en decodeer ze een voor een
                s = mod(r * checkMatrix',2); % bereken het syndroom
                for j=1:2^(n-k)
                    if(syndromes(j,:)==s) 
                        break; % zoek het syndroom op in de syndroomtabel
                    end
                end
                c = mod(r+cosets(j,:),2); % tel de gevonden cosetleider uit de syndroomtabel op bij het ontvangen codewoord om het oorspronkelijke codewoord te kunnen reconstueren
                bitdec((i-1)*k+1:i*k)=c(1:11); % het informatiewoord is bevat in de eerste 11 bits van het codewoord aangezien G in zijn systematische vorm staat
            end  
            
                        

        end

         % Functies voor de productcode
         
        function bitenc = Prod_encode(bitstring)
            % Functie die een bitstring encodeerd met de productcode.
            % input:
            % bitstring: vector van ongecodeerde bits met lengte 1xN
            
            bitstring = bitstring(:)';  % zorg ervoor dat de input een rijvector is
            N = length(bitstring);
            N_codewords = ceil(N/11/8);   
            
            % voeg nullen aan de bitstring toe indien N niet deelbaar is
            % door een geheel aantal codewoorden
            bitstring = [bitstring zeros(1, N_codewords*11*8-N)];
            bitenc = zeros(1, (8+1)*15*N_codewords);
            outputkolom = 1;
            % Product code hier
            for j = 1:N_codewords
                P_code_matrix = zeros(9,15);
                for i = 1:8
                    P_code_matrix(i,1:15) = Channel_Coding.Ham_encode(bitstring((j-1)*(8*11)+(i-1)*11+1:(j-1)*(8*11)+(i)*11));
                end
                for i = 1:15
                    P_code_matrix(9,i) = mod(sum(P_code_matrix(1:8,i)),2);
                end
                for i = 1:9
                    bitenc(1,outputkolom:outputkolom+14) = P_code_matrix(i,1:15);
                    outputkolom = outputkolom +15;
                end
            end

             
            
        end
        
        function bitdec = Prod_decode(bitenc)
            % Functie die een productcode decodeerd (foutcorrectie)
            % input:
            % bitenc: vector met gecodeerde bits: lengte moet deelbaar zijn door 15
            
            bitenc = bitenc(:)';    % zorg ervoor dat de input een rijvector is
            N = length(bitenc);
            N_codewords = N/15/(8+1);
            
            if(mod(N, 15*(8+1)) ~= 0)
                error('input is geen geheel aantal codewoorden.');
            end

            % output: de gedecodderde bits: lengte 11*8*N_codewords
            bitdec = zeros(1, 11*8*N_codewords);
            
            % errorcorrectie hier
            
            generator = [1 0 0 1 1]; 
            Hsys = Channel_Coding.Systematic_Check_Matrix(generator); 
            [syndromes,cosets]=Channel_Coding.Syndrome_Table(Hsys);
            
            bitenc = reshape(bitenc, 15,N_codewords*9)';
            N_Hamwords = N_codewords*8;
            Hamwords = zeros(N_Hamwords,15);
            
           % foutcorrectie
            for i=1:N_codewords     % i is het nummer van het productcodewoord
               % step 1
               matrix = bitenc((i-1)*9+1:i*9,:) ;    % productcodewoord
               parbits = mod(sum(matrix),2)  ;       % pariteitsbits voor alle kolommen van de matrix     
               syns = zeros(9,4);   % syndromen voor alle rijen van de matrix
               for j=1:9
                  syns(j,:)=mod(matrix(j,:)*Hsys',2); 
               end
               % step 2
               synsum = sum(syns')' ;    % som van elk syndroom in syns om te kijken of ze nul zijn
               if (sum(synsum)~=0)
                  % step 3
                  synindex = find(synsum);       % indexen van de syndromen in syns die niet nul zijn
                  
                  for j=1:length(synindex)  % j is de index voor de elementen in synindex => synindex(j)= index syndroom dat niet nul is in syns => syns(synindex(j),:)=niet-nul syndroom
                                            % => matrix(synindex(j),:) is foutief codewoord
      
                      foutpos = 0;      % positie van de fout in codewoord op rij synindex(j) => foutief bit =  matrix(synindex(j),foutpos)
                      for g=2:size(syndromes)
                          
                          if(syndromes(g,:)==syns(synindex(j),:))
                              foutpos = find(cosets(g,:));
                          end
                      end
                      
                      % step 4
                      if (parbits(foutpos)==1)
                          matrix(synindex(j),foutpos) = ~matrix(synindex(j),foutpos);
                          parbits(foutpos)=0;
                          syns(synindex(j),:) = [0 0 0 0];                         
                      end
 
                  end
                  
                  % einde step 4
                  
                  % step 5
                  synsum = sum(syns')'; % som van elk syndroom in syns om te kijken of er nog niet-nul syndromen over zijn
                  if (sum(synsum)~=0)
                      % er zijn nog niet-nul syndromen over => rijen of
                      % kolommen met meerdere fouten
                      synindex = find(synsum) ;      % indexen van de syndromen in syns die niet nul zijn
                      aantal = length(synindex) ;    % aantal niet-nul syndromen
                      if (aantal==1) 
                          matrix(synindex,:) = bitxor(matrix(synindex,:),parbits);
                          parbits = bitxor(parbits,parbits);
                          syns(synindex,:) = [0 0 0 0] ; 
                      elseif (aantal==2)
                          
                          if (syns(synindex(1),:)==syns(synindex(2),:))
                              % de 2 niet-nul syndromen zijn gelijk
                              foutpos = 0;      % positie van de fout in codewoord op rij synindex(1) => foutief bit =  matrix(synindex(1),foutpos)
                              for g=2:size(syndromes)

                                  if(syndromes(g,:)==syns(synindex(1),:))
                                      foutpos = find(cosets(g,:));
                                  end
                              end
                              
                              matrix(synindex(1),foutpos) = ~matrix(synindex(1),foutpos);
                              matrix(synindex(2),foutpos) = ~matrix(synindex(2),foutpos);
                              
                              syns(synindex(1),:) = [0 0 0 0];   
                              syns(synindex(2),:) = [0 0 0 0];   
                              
                          end
                          
                      end
                      
                  end
                  
                  
               end
                % einde correctie voor productcodewoord i
                Hamwords((i-1)*8+1:i*8,:) = matrix(1:8,:);
                              
            end
            % einde van de foutcorrectie 

            % decodering
            Hamwords = reshape(Hamwords',1, N_Hamwords*15);
            % output: de gedecodderde bits: lengte 11*8*N_codewords
            bitdec = Channel_Coding.Ham_decode(Hamwords);
                     
        end
        
        
    end
    
end