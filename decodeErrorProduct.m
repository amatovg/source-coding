%% 3.2.6

Homer

% Product code 
n=15;
k=11;
l=8;

% bit error probability
p =  [0.001 0.003 0.01 0.03 0.1 0.3];

% length of the bitsequence
length = 879995;
%length = 87995;
%length = 8795;
% number of informationbits
N_infowords= ceil(length/(k*l))
% bitsequence
bitstring = Channel_Coding.random_bitstring(length);

% add zeroes to the bitstring in case length isn't a multiple of a discrete
% number of codewords.
bitstring = [bitstring zeros(1, N_infowords*(k*l)-length)];  

% encode the bitsequence
bitencProd = Channel_Coding.Prod_encode(bitstring);
bitencHam = Channel_Coding.Ham_encode(bitstring);

% decode errors = number of informationwords (= l*k bits) that are different
N_decErrorsProd = zeros(1,size(p,2));
N_decErrorsHam = zeros(1,size(p,2));

% decode errors for the different error probabilities
for j=1:size(p,2)

    % introducing the same errors on both the encoded bitsequences
    sizeProd = size(bitencProd,2);
    sizeHam = size(bitencHam,2);
    errorsProd = (rand(1,sizeProd)<=p(j));
    errorsHam = errorsProd(1:sizeHam);
    biterrorProd = mod(errorsProd + bitencProd,2);
    biterrorHam = mod(errorsHam + bitencHam,2);

    % decoding the encoded bitsequence with errors
    bitdecProd = Channel_Coding.Prod_decode(biterrorProd);
    bitdecHam = Channel_Coding.Ham_decode(biterrorHam);

    % the difference between the decoded bitssequence and the original
    % bitsequence
    bitdiffProd = mod(bitdecProd+bitstring,2);
    bitdiffHam = mod(bitdecHam+bitstring,2);

    for i=1:(k*l):length
        if (sum(bitdiffProd(i:i+(k*l)-1))>0) 
            N_decErrorsProd(j) = N_decErrorsProd(j)+1;          
        end
        if (sum(bitdiffHam(i:i+(k*l)-1))>0) 
            N_decErrorsHam(j) = N_decErrorsHam(j)+1;          
        end
    end

end
N_decErrorsProd
N_decErrorsHam 

% decode error probability = decode errors / number of informationbits
Decode_Error_Probability_Prod = N_decErrorsProd/N_infowords
Decode_Error_Probability_Ham = N_decErrorsHam/N_infowords

% analytical decode error probability
Decode_Error_Analytical_Probability_Ham = 1-(1-1-((1-p).^15+15*p.*(1-p).^14)).^8

% plot of the simmulated decode error
% Hamming code: simulated and analytical, for l=8 codewords together
fig = figure;
loglog(p,Decode_Error_Probability_Ham,'-ob',p,Decode_Error_Analytical_Probability_Ham,':or','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63] );
grid on;
ylabel('Decode Error Probability');
xlabel('Bit Error Probability (p)');
title('Decode Error Probability for l=8 codewords of the (15,11) Hamming Code');
legend('Simulated','Analytical','Location','Southeast');
print(fig, '-djpeg', 'SimulatedDecodeErrorHamL.jpg');

% Product code and Hamming code simulated
fig = figure;
loglog(p,Decode_Error_Probability_Prod,'-ob',p,Decode_Error_Probability_Ham,':or','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63] );
grid on;
ylabel('Decode Error Probability');
xlabel('Bit Error Probability (p)');
title('Decode Error Probability for the Product Code and (15,11) Hamming Code');
legend('Product','Hamming','Location','Southeast');
print(fig, '-djpeg', 'SimulatedDecodeErrorProdHam.jpg');

% Product code simulated and Hamming code analytical
fig = figure;
loglog(p,Decode_Error_Probability_Prod,'-ob',p,Decode_Error_Analytical_Probability_Ham,':or','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63] );
grid on;
ylabel('Decode Error Probability');
xlabel('Bit Error Probability (p)');
title('Decode Error Probability for the simulated Product Code and analytical (15,11) Hamming Code');
legend('Product','Hamming','Location','Southeast');
print(fig, '-djpeg', 'SimulatedDecodeErrorProdSimHamAna.jpg');