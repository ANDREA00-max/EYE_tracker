% -----------------------------------------------------------------------
% SCRIPT PER IL RIDIMENSIONAMENTO DEI VETTORI ESTRATTI DAL MOBILE:
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

% -----------------------------------------------------------------------
% Caricamento file contenente gli istanti di inzio estrapolati dai video.
% -----------------------------------------------------------------------
load([pwd, '/Dati/istanti_inizio.mat']); %sto caricando un file che si chiama 
% istanti_inizio.mat che si trova dentro la cartella Dati però fuori dalle 
% cartelle dei soggetti. Attualmente non c'è, dovete inserirla voi
% Supponiamo che la variabile all'interno del file l'avete chiamato inizio
% ed è un vettore con 30 valori (uno per ogni soggetto considerato)

for i =1:length(path)
    
    % carico i dati del mobile
    load([path(i).folder '/' path(i).name '/' path(i).name '_mobile.mat'])
    
    % trovo il campione a cui corrisponde l'inizio della baseline
    istante = find(Time_M > istanti_inizio(i),1,'first');
    
    % elimino tutti i segmenti inziali di tutti i vettori
    BEventInfo_M = BEventInfo_M(istante : end);
    LPupilDiametermm_M = LPupilDiametermm_M(istante : end);
    RPupilDiametermm_M = RPupilDiametermm_M(istante : end);
    Time_M = Time_M(istante : end);
    
    
    % salvo in un altro file all'interno della stessa cartella le variabili
    % aggiornate
    save([path(i).folder '/' path(i).name '/' path(i).name '_mobile_aggiornato.mat'],...
        'Time_M', 'LPupilDiametermm_M', 'RPupilDiametermm_M', 'BEventInfo_M')
    
end
    


