path=dir([pwd,'/Dati/S*']);

for i=1:length(path)
    % carico cartella con tutti i dati di interesse del soggetto,12:
    % RM,PSD_RM,LM,PSD_LM;RF,PSD_RF,LF,PSD_LF,// tacoframma, PSD_TG,
    % respirogramma, PSD_RG.
    %bisogna scelgiere se inserire i dati già filtrati e ricampionati o in
    %orgiinale dopo la ricostruzione.
    %dentro a plot c'e il nimero del soggetto e la posizione nella matrice
    %figure che viene creata da subplot, 
load([path_dati(i).folder '/' path_dati(i).name '/' path_dati(i).name 's'])

for j = 1:8
    subplot(6,2,1);
    plot(LPupilDiametermm_F_Final);
    
    subplot(6,2,2);
    plot(PSDLF);
    
    subplot(6,2,3);
    plot(LPupilDiametermm_M_Final);
    
    subplot(6,2,4);
    plot(PSD-LM);
    
    subplot(6,2,5);
    plot(RPupilDiametermm_F_Final);
    
    subplot(6,2,6);
    plot(PSD-RF);
    
    subplot(6,2,7);
    plot(RPupilDiametermm_M_Final);
    
    subplot(6,2,8);
    plot(PSD-RM);
    
    subplot(6,2,9);
    plot(sig_respirogramma);
    
    subplot(6,2,10);
    plot(PSD_RESP);
    
    subplot(6,2,11);
    plot(sig_tacogramma);
    
    subplot(6,2,12);
    plot(PSD_taco);
    
    
end

end