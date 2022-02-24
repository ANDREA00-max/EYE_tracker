% -----------------------------------------------------------------------
% SCRIPT PER LA SINCRONIZZAZIONE DEI VETTORI ESTRATTI DAL MOBILE E DAL FISSO:
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
     
    %CARICO I DATI DEL FISSO 
    load([path(i).folder '/' path(i).name '/' path(i).name '_fisso_ripulito.mat'])
    LPupilDiametermm_F_Finale=LPupilDiametermm_F;
    s1=medfilt1(LPupilDiametermm_F(10000:end-10000));
    RPupilDiametermm_F_Finale=RPupilDiametermm_F;
    
    %CARICO I DATI DEL MOBILE 
    load([path(i).folder '/' path(i).name '/' path(i).name '_mobile_aggiornato_ripulito.mat'])
    LPupilDiametermm_M_Finale=LPupilDiametermm_M;
    s2=	medfilt1(LPupilDiametermm_M(10000:end-10000));
    RPupilDiametermm_M_Finale=RPupilDiametermm_M;
    
    %UTILIZZO LA CROSS-CORRELAZIONE PER OTTENERE IL DELAY TEMPORALE
    
    [C21,lag21] = xcorr(s1,s2);
    C21 = C21/max(C21);
    
    
    [M21,I21] = max(C21);
    t21 = lag21(I21);
    
    

    
    %RINOMINO I RISULTATI
    
    LPupilDiametermm_F_Finale=LPupilDiametermm_F_Finale(t21:end);
    RPupilDiametermm_F_Finale=RPupilDiametermm_F_Finale(t21:end);
  
    
   %aggiorno i vettori temporali: lunghi come quelli del diametro e poi
   %fratto 30 perchè acquisiamo a 30 Hz
    tempo_F = (1:length(RPupilDiametermm_F_Finale)) / 30;  
    tempo_M = (1:length(RPupilDiametermm_M_Finale)) / 30;
        
    save([path(i).folder '/' path(i).name '/' path(i).name '_sincronizzati_xcorr'],...
         'LPupilDiametermm_F_Finale','RPupilDiametermm_F_Finale','RPupilDiametermm_M_Finale','LPupilDiametermm_M_Finale','tempo_M','tempo_F')

     figure
     plot(RPupilDiametermm_F_Finale)
     hold on
     plot(RPupilDiametermm_M_Finale)
     ylim([1 6])


end