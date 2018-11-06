function err = compute_metrics(Y,Y_hat)

% COMPUTE_METRICS returns cell array of metrics between GT SRI and estimation
% err = COMPUTE_METRICS(Y,Yhat) returns cell array err of metrics btw Y and Y_hat
% INPUT ARGUMENTS:
%     Y: GT SRI
%     Y_hat: estimated SRI
% OUTPUT ARGUMENTS:
%     err: contains R-SNR, CC, SAM, ERGAS

err = {r_snr(Y,Y_hat), cc(Y_hat), sam(Y,Y_hat), ergas(Y,Y_hat,1/4)};

end

