%% nuovo calcolo HF con integrale adattivo

%sfrutto la funzione findpeaks che mi trova i max relativi della funzione e
%i suoi indici, nelle LF integrero sul primo max relativo , nelle HF sul
%secondo max relativo, sarebbe da implementare sui soggetti in cui il picco
%delle LF è nell'elemento 1 e non viene trovato dalla funzione come max
%relativo

path=dir([pwd,'/Dati/S*']);



MAT_LOW= [];
MAT_HI= [];
MAT_RATIO = [];
MAX_vec= [];
MAX_vec_2= [];
MAT_SUM=[];
MIN_vec=[];
for j=1:length(path)
     
    %CARICO I DATI DEL FISSO 
row_low = [];
row_hi = [];
row_ratio= [];
row_sum = [];
   load([path(j).folder '/' path(j).name])
for i = 1:9
% MAX = find (statsLM(i).cPSDtotTimeFrames{1, 1}==...
%     max(statsLM(i).cPSDtotTimeFrames{1, 1}))
% 
% MAX_vec = [MAX_vec;MAX];
% MAX_2 = find (statsLM(i).cPSDtotTimeFrames{1, 1}(30:end)==...
%     max(statsLM(i).cPSDtotTimeFrames{1, 1}(30:end)))
% 
% MAX_vec_2 = [MAX_vec_2;MAX_2];
% 
% MIN =find (statsLM(i).cPSDtotTimeFrames{1, 1}(1:50)==...
%     min(statsLM(i).cPSDtotTimeFrames{1, 1}(1:50)))
% MIN_vec = [MIN_vec;MIN];
[pks,locs] = findpeaks(statsLM(i).cPSDtotTimeFrames{1, 1})




  xLF = statsLM(i).PSDfreqs(1:locs(1)+10);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (1:locs(1)+10);

 ILF = trapz(xLF,yLF);
 

 
 xHF = statsLM(i).PSDfreqs(locs(2)-10:locs(2)+10);
 yHF = statsLM(i).cPSDtotTimeFrames{1,1} (locs(2)-10:locs(2)+10);

 IHF = trapz(xHF,yHF);
 

 
 row_low = [row_low,ILF];
 row_hi = [row_hi,IHF];
 row_ratio= [row_ratio, ILF/IHF];
 row_sum= [row_ratio, ILF+IHF];
end
MAT_LOW= [MAT_LOW;row_low];
MAT_HI= [MAT_HI;row_hi];
MAT_RATIO = [MAT_RATIO;row_ratio];
MAT_SUM = [MAT_SUM;row_sum];
end

