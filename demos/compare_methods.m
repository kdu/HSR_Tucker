function [res] = compare_methods(SRI, HSI, MSI, DegMat, ds, methods)
  res = cell(size(methods,1), 6);
  for i=1:size(methods,1)
    % Initialize parameters and find name
    opt = struct();
    blockstr = '';
    if size(methods,2) > 3 && ~isempty(methods{i,4})
      opt.Nblocks = methods{i,4};
      blockstr = sprintf('(%d,%d) ', opt.Nblocks(1), opt.Nblocks(2)); 
    end  
    methname = sprintf('%s %s%s', methods{i,1}, blockstr, methods{i,3});
    fprintf('Running method %s\n', methname);
    
    Params = sort(fieldnames(DegMat));  
    ParamStr = sprintf('DegMat.%s, ', Params{:});
    tic,
    eval(sprintf('[Y_hat,~] = %s(HSI, MSI, %s %s, opt);', ...
             methods{i,2}, ParamStr,  methods{i,3}));
    time = toc;
    res(i,:) = [methname, compute_metrics(SRI,Y_hat,ds(1),ds(2)), time];
  end       
end



