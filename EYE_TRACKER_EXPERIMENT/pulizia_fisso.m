%%pulizia_segnale_fisso se esiste un istante in cui il software dell'eye  
%tracker ha riconosciuto un blink, setto il valore corrispondente a quell'istante come NaN

close all;
clear;
clc;


path=dir([pwd,'/Dati/S*']);

for i=1:length(path)
     
    %%carico un soggetto diverso ad ogni iterazione del ciclo for
    
    load([path(i).folder '/' path(i).name '/' path(i).name '_fisso.mat'])
  
     %%"LBlink" e "RBlink" sonp delle variabili logiche che mi indicano dove sono i Blink nel vettore dei dati,
    %e sono date dalla variabile categorica "LEventInfo_F","REeventInfo_F".
    
    L_Blink=LEventInfo_F=='Blink';
    R_Blink=REventInfo_F=='Blink';
    
    
    %Pongo uguali a NaN i valori agli istanti riconosciuti dalla variabile
    %Blink. Elimino inoltre i valori non fisiologici (D>7, D<2) che sono
    %dovuti ad errori di lettura dell'Eye Tracker.
    
    
    LPupilDiametermm_F(L_Blink)=NaN;
    RPupilDiametermm_F(R_Blink)=NaN;
    LPupilDiametermm_F(LPupilDiametermm_F<1 | LPupilDiametermm_F>7)=NaN;
    RPupilDiametermm_F(RPupilDiametermm_F<1 | RPupilDiametermm_F>7)=NaN;
    
    %Eseguo il resampling per eliminare l'errore di acquisizione per il
    %quale la frequenza di campionamento non era esattamente costante.
    %'pchip' ricostruisce i tratti di segnale "danneggati" dai NaN. NB: i
    %primi istanti danno dei problemi con i valori che schizzano altissimi
    %quindi per visualizzare il grafico corretto settate che la y vada ad
    %esempio da 0 a 6mm, altrimenti non si vede nulla
    
    
    [L_Resample,t]=resample(LPupilDiametermm_F, Time_F, 30, 'pchip');
    [R_Resample,t]=resample(RPupilDiametermm_F, Time_F, 30, 'pchip');
    
    
    
    LPupilDiametermm_F=L_Resample;
    RPupilDiametermm_F=R_Resample;
    
    Time_F = t;
     %Salvo i risultati ottenuti in un nuovo file Mat.
     
     
    save([path(i).folder '/' path(i).name '/' path(i).name '_fisso_ripulito.mat'],...
        'Time_F', 'LPupilDiametermm_F', 'RPupilDiametermm_F', 'LEventInfo_F', 'REventInfo_F')
end