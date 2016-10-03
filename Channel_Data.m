function post_data = Channel_Data(data,p)
% Function that simulates a noisy, binary, symmetric channel which is
% capable of inducing errors in the data who pass through.
% input:
% data: a bit sequence at the beginning of the channel
% p: error probability: chance for each bit of getting complemented 
%    (inducing an error)
% output:
% post_data: the data at the end of the channel (modified)

length = size(data,2);
errors = (rand(1,length)<=p);
post_data = mod(errors + data,2);

end