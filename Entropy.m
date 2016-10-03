%% 3.1.3 Entropy

% load 'Xfile.mat'
load('Xfile2.mat');


% number of bits in 1 macrosymbol is K
% number of macrosymbols is 2^K
% H(x) = - sum( rel_freq(x) * log2( rel_freq(x) ) )   [bits/macrosymbol]
% E(n) = 1/K sum ( rel_freq(x) * n(x) )     [codebits/source symbol]


% entropy for K=1..10
entropy = zeros(1,10);      % for Xfile
entropy_uniform = zeros(1,10);      % for uniform frequencies

% average number of codebits per source symbol
meanLength = zeros(1,10);

% max recursion needs to be high enough
set(0,'RecursionLimit',1500)

for K=1:10
    % counter holds the number of times a macrosymbol is in the image
    counter=zeros(1,2^K);
    
    % x is a macrosymbols from the data
    x = zeros(1,K);

    % counting the macrosymbols
    for i=1:K:length(bitsrc)
        x = bitsrc(i:i+K-1);
        index = bi2de(x);      
        counter(index+1) = counter(index+1)+1;     
    end
    
    % relative frequencies of the macrosymbols
    rel_freq=counter/sum(counter);
    rel_freq_uniform = ones(1,2^K)/(2^K);
    
    % calculate entropy for K
    mask=rel_freq(rel_freq>0); % to avoid NaN if a macrosymbol doesn't occur
    entropy(K)= - sum(mask.*(log2(mask)));
    entropy_uniform(K) = -sum(rel_freq_uniform.*(log2(rel_freq_uniform)));
    
    % Huffman
    alphabet_K = cell(1,2^K);
    codewords_K = Source_Coding.create_codebook(alphabet_K, rel_freq);
    lengthCodewords = cellfun(@length,codewords_K);
    meanLength(K) = sum(lengthCodewords.*rel_freq)/K;
    
end

% plot of the entropy
x1 = [1:10];    % values of K
fig = figure;
plot(x1,entropy,'b',x1,entropy_uniform,'r');
ylabel('Entropy');
xlabel('Macrosymbol size K (# bits)');
title('Entropy for macrosymbols');
legend('not uniform','uniform','Location','Southeast');
print(fig, '-djpeg', 'Entropy.jpg');


%% 3.1.4 E(n)

% E(n) = 1/K sum ( rel_freq(x) * n(x) )     [codebits/source symbol]
% limits:  H(x) <= E(n) < H(x)+1/K   with everything in [codebits/source symbol] !!!!

% entropy per source symbol 
entropySource = entropy./x1;

% boundaries
upperBound = entropySource+1./x1;
lowerBound = entropySource;

% plot of the average lenght
fig = figure;
plot(x1,upperBound,'r',x1,meanLength,'g',x1,lowerBound,'b')
ylabel('Mean Length E[n]')
xlabel('Macrosymbol size K (#bits)');
title('Mean Length for Source Symbols');
legend('Upper Bound','Mean Length','Lower Bound','Location','Northeast');
print(fig, '-djpeg', 'MeanLenght.jpg');






%expected value E[n] for K=1..10


