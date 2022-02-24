
path=dir([pwd,'/Dati/S*']);

for i=1:length(path)
load([path(i).folder '/' path(i).name '/' path(i).name '_mobile_aggiornato.mat'])
k=0; 
% calcolo il differnziale
LDIFF = abs(diff(LPupilDiametermm_M));
RDIFF = abs(diff(RPupilDiametermm_M));
 LPupilDiametermm_M_DIFF=LPupilDiametermm_M;
 RPupilDiametermm_M_DIFF=RPupilDiametermm_M;
% pulizia mobile 
% il diff avrà un elemento in meno, scarto diff>0.375
% 
for k = 3:(length(Time_M)-5)
    if (LDIFF(k)>0.375 )
       LPupilDiametermm_M_DIFF(k-2:k+4)=NaN;
    
    end
    
    if  LPupilDiametermm_M_DIFF(k)==0
        LPupilDiametermm_M_DIFF(k-3:k+3)=NaN;
    end

end



for k = 1:(length(Time_M)-1)
    if (RDIFF(k)>0.375 || LDIFF(k)==0) 
       RPupilDiametermm_M_DIFF(k+1)=NaN;
    end


    if  LPupilDiametermm_M_DIFF(k)==0
        LPupilDiametermm_M_DIFF(k-3:k+3)=NaN;
    end
end

 save([path(i).folder '/' path(i).name '/' path(i).name 'mobile_aggiornato_DIFF'],...
 'RPupilDiametermm_M_DIFF', 'LPupilDiametermm_M_DIFF' ,'Time_M', 'BEventInfo_M')
end

%% nello script pulizia mobile bisogna aggiungere_DIFF alle varibili
