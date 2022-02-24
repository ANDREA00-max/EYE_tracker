% -----------------------------------------------------------------------
% SCRIPT PER IL SALVATAGGIO DEI RISULTATI OTTENUTI CON LE DIVERSE FUNZIONI
% DI SINCRONIZZAZIONE
% 
%
%

close all;
clear;
clc;

% -----------------------------------------------------------------------
% Inizializzazione variabli.
% -----------------------------------------------------------------------

path=dir([pwd,'/Dati/S*']);

for i=1:length(path)
     
    %CARICO I DATI SINCRONIZZATI CON XCORR e con TIMEDELAY, per s03 ed s12
    %salvo i valori ottenuti con timedelay, altrimenti salvo quelli
    %ottenuti con xcorr
    
    if (i==2 | i==11)
        load([path(i).folder '/' path(i).name '/' path(i).name '_sincronizzati.mat'])
        LPupilDiametermm_F_Final=LPupilDiametermm_F_Finale;
        LPupilDiametermm_M_Final=LPupilDiametermm_M_Finale;
        RPupilDiametermm_F_Final=RPupilDiametermm_F_Finale;
        RPupilDiametermm_M_Final=RPupilDiametermm_M_Finale;
        tempo_Fisso=tempo_F;
        tempo_Mobile=tempo_M;
    else
    load([path(i).folder '/' path(i).name '/' path(i).name '_sincronizzati_xcorr.mat'])
        LPupilDiametermm_F_Final=LPupilDiametermm_F_Finale;
        LPupilDiametermm_M_Final=LPupilDiametermm_M_Finale;
        RPupilDiametermm_F_Final=RPupilDiametermm_F_Finale;
        RPupilDiametermm_M_Final=RPupilDiametermm_M_Finale;
        tempo_Fisso=tempo_F;
        tempo_Mobile=tempo_M; 
    end
    save([path(i).folder '/' path(i).name '/' path(i).name '_sincronizzazione_finale'],...
         'LPupilDiametermm_F_Final','RPupilDiametermm_F_Final','RPupilDiametermm_M_Final','LPupilDiametermm_M_Final','tempo_Mobile','tempo_Fisso')
end
        