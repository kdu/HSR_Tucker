function [SRI_hat,info] = hysure_adaptor(HSI, MSI, Pm, R, opts)
%Changes the order of parameters 
  [SRI_hat,info] = hysure(HSI,MSI, opts);
end

