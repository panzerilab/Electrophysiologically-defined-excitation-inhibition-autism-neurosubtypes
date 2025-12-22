function len_pool = automatic_pooling(data,fs)
%   function to detect the length for pooling
%   check for the the lag at which the autocorrelation of the signal drops
%   below significant level, if that happens less than 1 minute, it returns
%   1 minute 

n_ch = size(data,2);
T_1m = 60*fs;
n_lags = round(2*60*fs);
pooling_len = zeros(n_ch,1);
for ch=1:n_ch
    [acf,~,bounds] = autocorr(data(:,ch),n_lags);
    pooling_len(ch) = find(abs(acf)<bounds(1),1);
end

tmp_pool = max(pooling_len);
len_pool = max(tmp_pool,T_1m);

while rem(size(data,1),len_pool) ~= 0
      len_pool = len_pool + 1;
end

end

