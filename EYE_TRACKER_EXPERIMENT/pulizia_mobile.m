%pulizia segnale mobile: se esiste un istante in cui il software dell'eye
%tracker ha riconosciuto un blink, setto il valore corrispondente a quell'istante come NaN

close all;
clear;
clc;


path=dir([pwd,'/Dati/S*']);

for i=1:length(path)
     
    %%carico un soggetto diverso ad ogni iterazione del ciclo for
    
    load([path(i).folder '/' path(i).name '/' path(i).name '_mobile_aggiornato.mat'])
    
    %%"Blink" è una variabile logica che mi indica dove sono i Blink nel vettore dei dati,
    %ed è data dalla variabile categorica "BEventInfo_M".
    
    Blink=BEventInfo_M=='Blink';
    
    %Pongo uguali a NaN i valori agli istanti riconosciuti dalla variabile
    %Blink. Elimino inoltre i valori non fisiologici (D>7, D<2) che sono
    %dovuti ad errori di lettura dell'Eye Tracker.
    
    LPupilDiametermm_M(Blink)=NaN;
    LPupilDiametermm_M(LPupilDiametermm_M<1 | LPupilDiametermm_M>7)=NaN;
    RPupilDiametermm_M(Blink)=NaN;
    RPupilDiametermm_M(RPupilDiametermm_M<1 | RPupilDiametermm_M>7)=NaN;
    
    %Eseguo il resampling per eliminare l'errore di acquisizione per il
    %quale la frequenza di campionamento non era esattamente costante.
    %'pchip' ricostruisce i tratti di segnale "danneggati" dai NaN. NB: i
    %primi istanti danno dei problemi con i valori che schizzano altissimi
    %quindi per visualizzare il grafico corretto settate che la y vada ad
    %esempio da 0 a 6mm, altrimenti non si vede nulla
    
    [L_ricampionato,t]=resample(LPupilDiametermm_M,Time_M,30,'pchip');
    
    LPupilDiametermm_M=L_ricampionato;
    
    [R_ricampionato,t]=resample(RPupilDiametermm_M,Time_M,30,'pchip');
    
    RPupilDiametermm_M=R_ricampionato;
    
    Time_M = t;
    
    %Salvo i risultati ottenuti in un nuovo file Mat.
    
    save([path(i).folder '/' path(i).name '/' path(i).name '_mobile_aggiornato_ripulito.mat'],...
        'Time_M', 'LPupilDiametermm_M', 'RPupilDiametermm_M', 'BEventInfo_M')
end