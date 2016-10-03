classdef Source_Coding
    
    methods(Static=true)
        
        % Functies voor de Huffmancode op te stellen
        function codewords = create_codebook(alphabet, rel_freq)
            % Functie die de Huffmancode opstelt voor gegeven alfabet en
            % relative frequenties
            % input:
            % alphabet: 1xN cell array vb alphabet = {'A', 'B', 'C', 'D'}
            % rel_freq: relative frequenties voor elke letter in het
            % alfabet: 1xN vector vb: rel_freq = [0.5 0.1 0.3 0.1];
            
            N = numel(alphabet);
            rel_freq = rel_freq/sum(rel_freq);
            
            %output the codewords
            % example: codewords = {[0], [1 1 0], [1 0], [1 1 1]};
            codewords = cell(1, N);
            
            % afronden van de frequenties voor zelfde code als de handmatige code te
            % bekomen
            rel_freq= round(rel_freq.*10000);
            
            % kolom1 = relatieve frequenties, kolom2 -index van de
            % macrosymbolen (negatief voor aflopende sortering te bekomen voor de index)
            matrix = [rel_freq' -1.*[1:N]'];
            
            % binaire boom voor de huffmancode in matrixvorm: kolom1 =
            % som van de frequenties van de macrosymbolen van de deelboom,
            % kolom2-17 = indexen van de macrosymbolen op deeindtakken van de deelboom
            % boom is gesorteerd van kleinste frequentie naar hoogste
            % frequentie en dan van hoogste index naar laagste index zodat
            % bovenste rij bit 1 krijgt en tweede rij bit 0 (omgekeerde
            % volgorde dan bij de handmatige boom)
            % N+2 kolommen: 1 kolom frequentie, N kolommen nodig voor
            % indexen op het einde van het algoritme, 1 kolom met 0 voor
            % het einde van de indexen te vinden
            tree = zeros(N,N+2);
            tree(:,1:2) = sortrows(matrix);
            tree = abs(tree);
            
            % algoritme:
            % bit 1 of 0 toevoegen aan de codewoorden met index in de eerste twee rijen,
            % indexen van de eerste rij toevoegen aan de tweede rij,
            % frequentie van de eerste rij optellen bij de frequentie van de tweede rij
            % eerste rij verwijderen
            % boom opnieuw sorteren
            % herhalen tot er maar 1 rij overblijft
            while(size(tree,1)>1)
                % bit 1 voor eerste rij
                % bit 0 voor tweede rij
                min0 = tree(2,1);
                min1 = tree(1,1);
                
                einde0 = find(tree(2,:)==0,1);
                for i=2:einde0-1
                    index0 = tree(2,i);
                    codewords(index0) = strcat('0',codewords(index0));
                end
                
                einde1 = find(tree(1,:)==0,1);
                for j=2:einde1-1
                    index1 = tree(1,j);
                    codewords(index1) = strcat('1',codewords(index1));
                    tree(2,einde0+j-2) = index1;
                end
                
                tree(2,1) = min0+min1;
                tree = tree(2:end,:);
                tree(:,2:end) = -1.*tree(:,2:end);
                tree = sortrows(tree);
                tree = abs(tree);
                
            end
            for i=1:N
                temp = codewords{i};
                codewords{i}=str2num(temp(:))';
            end
            
        end
        
        function codewords = create_canonical_codebook(alphabet, lengths)
            % Functie die de Cononical Huffmancode opstelt voor gegeven alfabet en
            % codelengtes
            % input:
            % alphabet: 1xN cell array vb alphabet = {'A', 'B', 'C', 'D'}
            % lengths: lengtes van elk codewoord: [1 3 2 3];
            
            N = numel(alphabet);
            % the first macrosymbol (i.e. the first letter of the alpabet) receives
            % index 1, the second macrosymbol receives index 2,...
            indexAlphabet = [1:N];
            
            % matrix with in the first column the indices of the macrosybols in
            % alphabetical order, in the second column the lengths of the codewords,
            % and in the last column the not yet defined codewords
            % matrixSorted is the same matrix but first sorted by length of the
            % codewords and then sorted alphabetically
            matrix = [indexAlphabet' lengths' zeros(N,1)];
            matrixSorted = sortrows(matrix,[2,1]);
            
            % algorithm for the canonical huffmancode with the codewords in decimal
            for i=2:N
                matrixSorted(i,3) = matrixSorted(i-1,3)+1;
                x = de2bi(matrixSorted(i,3));
                while (length(x)<matrixSorted(i,2))
                    matrixSorted(i,3)=matrixSorted(i,3)*2;
                    x = de2bi(matrixSorted(i,3));
                end
            end
            
            matrixSorted = sortrows(matrixSorted,1);
            
            % output the codewords
            % example: codewords = {[0], [1 1 0], [1 0], [1 1 1]};
            codewords = cell(1, N);
            % the smallest codeword a number of zeros with the number
            % defined by the length of the codeword
            codewords{1} = zeros(1,matrixSorted(1,2));
            for i=2:N
                codewords{i} = de2bi(matrixSorted(i,3),'left-msb');
            end
            
        end
        
        
        
        % Functies voor het encoderen/decoderen
        function compressed = Huffman_encode(data, alphabet, codewords)
            % Functie die een bitstring comprimeerd met een Huffmancode
            % input:
            % data: de data die gecomprimeerd moet worden
            % alphabet: cell array of rijvector met alle mogelijke symbolen
            % codewords: cell array met de codewoorden voor elke letter in
            % het alfabet.
            
            N = length(data);
            N_symbols = length(alphabet);
            
            if iscell(alphabet) && max(cellfun(@length, alphabet))>1
                error('ERROR: een symbool uit het alfabet mag maar 1 cijfer of alfanumerieke waarde bevatten');
            end
            
            if ischar(data)
                % converteer de alfanumerieke waardes naar numeriek waardes
                numeric_data = double(data);
                numeric_symbols =  cellfun(@double, alphabet);
            elseif isrow(data)
                numeric_data = data;
                if iscell(alphabet),    numeric_symbols = cell2mat(alphabet);
                else numeric_symbols = alphabet;
                end
            else
                error('ERROR: input moet ofwel een rijvector of een karakater-string zijn');
            end
            
            % De kern van deze functie: vervang elke letter uit het alfabet
            % met het corresponderend codewoord
            compressed = cell(1,N);
            for n=1:N_symbols
                idx = numeric_data==numeric_symbols(n);
                compressed(idx) = {codewords{n}};
            end
            
            % Als er nog lege cellen overblijven wil dit zeggen dat het
            % woordenboek niet alle mogelijke symbolen uit de data heeft.
            if nnz(cellfun(@isempty,compressed))
                error('ERROR: woordenboek incompleet');
            end
            
            % zet om naar een rijvector
            compressed = cell2mat(compressed);
            
        end
        
        function decompressed = Huffman_decode(data_compressed, alphabet, codebook)
            % Functie die een gecomprimeerde bitsequentie decomprimeerd
            % Optimaal werkt deze functie met een boom maar hier
            % decomprimeren we simpelweg door te itereren.
            % input
            % data_compressed
            % alphabet: cell array of rijvector met alle mogelijke symbolen
            % codewords: cell array met de codewoorden voor elke letter in
            % het alfabet.
            
            N = length(data_compressed);
            
            % voorbereiding: splits de codewoorden op per lengte
            lengths = cellfun(@length, codebook);
            max_len = max(lengths);
            
            cl = cell(1,max_len);
            idxl = cell(1,max_len);
            for l=1:max_len
                idx_l = find(lengths==l);
                % zet de codewoorden van lengte l om naar hun decimale vorm
                cl{l} = cellfun(@bi2de, codebook(idx_l));
                idxl{l} = idx_l;
            end
            
            % decompressie: doorloop de data en zoek naar codewoorden
            output = [];    % bevat de indices van de letters die herkend zijn.
            idx = 1;        % index in de data
            window = 1;     % bereik waarover we zoeken: we beginnen met zoeken naar codewoorden van lengte 1
            while idx+window <= N+1
                idx_code = find(cl{window} == bi2de(data_compressed(idx:idx+window-1)));
                if ~isempty(idx_code)
                    %codewoord gevonden, voeg dit toe aan de output en ga verder
                    output = [output idxl{window}(idx_code)];
                    idx = idx + window;
                    window = 1;
                else
                    % nog geen codewoord gevonden; vergroot het bereik
                    window = window+1;
                end
            end
            
            % zet de output om naar letters uit het alfabet
            decompressed = alphabet(output);
            
        end
    end
    
end