%% 3.2.3
% (15,11) Hamming code
n=15;
k=11;

% bit error probability
p =  [0.001 0.003 0.01 0.03 0.1 0.3];

% length of the bitsequence
length = 109995;

% number of informationbits
N_infowords = ceil(length/k)
% bitsequence
bitstring = Channel_Coding.random_bitstring(length);

% encode the bitsequence
bitenc = Channel_Coding.Ham_encode(bitstring);

% decode errors = number of informationwords (= k bits) that are different
N_decErrors = zeros(1,size(p,2));

% add zeroes to the bitstring in case length isn't a multiple of a discrete
% number of codewords.
bitstring = [bitstring zeros(1, N_infowords*k-length)];  


% decode errors for the different error probabilities
for j=1:size(p,2)

    % introducing the errors on the encoded bitsequence
    biterror = Channel_Data(bitenc,p(j));

    % decoding the encoded bitsequence with errors
    bitdec = Channel_Coding.Ham_decode(biterror);

    % the difference between the decoded bitssequence and the original
    % bitsequence
    bitdiff = mod(bitdec+bitstring,2);

    for i=1:k:length
        if (sum(bitdiff(i:i+k-1))>0) 
            N_decErrors(j) = N_decErrors(j)+1;
        end
    end

end

N_decErrors

% decode error probability = decode errors / number of informationbits
DecodeErrorProbability = N_decErrors/N_infowords

% analytical decode error probability
% approximationProbability = 105*p.^2
DecodeErrorAnalyticalProbability = 1-((1-p).^15+15*p.*(1-p).^14)

% plot of the simmulated decode error
fig = figure;
loglog(p,DecodeErrorProbability,'-ob',p,DecodeErrorAnalyticalProbability,':or','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63] );
grid on;
ylabel('Decode Error Probability');
xlabel('Bit Error Probability (p)');
title('Decode Error Probability for (15,11) Hamming Code');
legend('Simulated','Analytical','Location','Southeast');
print(fig, '-djpeg', 'SimulatedDecodeErrorHamming2.jpg');

