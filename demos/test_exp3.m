% EXPERIMENT 3: COMPARISON OF METRICS WITH STEREO %
% Copyright (c) 2018 Clemence Prevost, Konstantin Usevich, Pierre Comon, David Brie
% https://github.com/cprevost4/HSR_Tucker
% Contact: clemence.prevost@univ-lorraine.fr


% %Salinas
% SRI = cell2mat(struct2cell(load('Salinas.mat')));
% SRI = crop(SRI,[80,84,size(SRI,3)]);
% SRI(:,:,[108:112 154:167 224]) = []; %Regions of water absorption (Salinas)

%Indian Pines
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption (Indian pines)
SRI(1,:,:) = []; SRI(:,1,:) = [];

Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

methods_nb = {'stereo','scott'};
names_nb = {'STEREO F=20','STEREO F=30','STEREO F=50','STEREO F=100','NaN',...
    'SCOTT = [14,14,15]','SCOTT = [24,24,25]',...
    'SCOTT = [30,30,16]','SCOTT = [40,40,6]','SCOTT = [58,58,6]'};
ranks_nb = {'20','[14,14,15]';
            '30','[24,24,25]';
            '50','[30,30,16]';
            '100','[40,40,6]';
            '0','[58,58,6]';
            };       
c = 0;
        
 for i=1:size(methods_nb,2)
     for j=1:size(ranks_nb,1)
         c = c+1;
         if ranks_nb{j,i}~="0"
            tic,
            eval(sprintf('[Y_hat,info] = %s(HSI, MSI, P1, P2, Pm, %s);', methods_nb{i},ranks_nb{j,i}));
            time = toc;
            eval('err = [compute_metrics(SRI,Y_hat,d1,d2), time];')
            res{1,c} = cell2mat(err)';
         else
            res{1,c} = NaN*ones(1,5)';
         end
     end
 end       
 save_text_table('res_exp3.txt', names_nb, cell2mat(res));



