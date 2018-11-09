% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr

addpath ../utils
addpath ../metrics
addpath ../methods
clear res

% Indian Pines
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
SRI(1,:,:) = []; SRI(:,1,:) = [];
% 
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
 
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

 methods_b = {'bscott'};
 names_b = {'BSCOTT R=[40,40,6]','BSCOTT R=[60,60,6]', 'BSCOTT R=[100,100,6]'};
 ranks_b = {'[40,40,6]';'[60,60,6]';'[100,100,6]'};
 c = 0;

 for i=1:size(methods_b,2)
     for j=1:size(ranks_b,1)
         c = c+1;
         if ranks_b{j,i} == "[100,100,4]"
             opts.Nblocks = [1,1];
         end
         if ranks_b{j,i}~="0"
            tic,
            eval(sprintf('[Y_hat,info] = %s(MSI, HSI, Pm, %s, opts);', methods_b{i},ranks_b{j,i}));
            time = toc;
            eval('err = [compute_metrics(SRI,Y_hat,d1,d2), time];')
            res{1,c} = cell2mat(err)';
         else
            res{1,c} = NaN*ones(1,5)';
         end
     end
 end  
 save_text_table('res_exp3_ku.txt', names_b, cell2mat(res));

