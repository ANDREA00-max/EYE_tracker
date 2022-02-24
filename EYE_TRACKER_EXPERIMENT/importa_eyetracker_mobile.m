% -----------------------------------------------------------------------
% SCRIPT PER L'IMPORTAZIONE DEI DATI:
% dalla cartella contenente i file di testo, lo script converte i dati in 
% variabili MATLAB e le salva.
% 
%
%
close all;
clear;
clc;

% -----------------------------------------------------------------------
% Inizializzazione variabli.
% -----------------------------------------------------------------------
folder = uigetdir;
path=dir([folder,'/*.txt']);
delimiter = '\t'; 
startRow = 34;

% -----------------------------------------------------------------------
% Formato per ogni linea del testo:
% -----------------------------------------------------------------------

formatSpec = '%f%C%f%f%f%f%f%f%f%f%f%C%[^\n\r]';

% -----------------------------------------------------------------------
% Apertura dei file in ciclo FOR.
% -----------------------------------------------------------------------

for i=1:length(path)
    
    filename = [path(i).folder '\' path(i).name];
    fileID = fopen(filename,'r');

    %%% Lettura dei dati.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    %%% Chiusura del file.
    fclose(fileID);

    %%% Allocate imported array to column variable names
    Time_M = dataArray{:, 1};
    Type_M = dataArray{:, 2};
    Trial_M = dataArray{:, 3};
    LDiaXpx_M = dataArray{:, 4};
    LDiaYpx_M = dataArray{:, 5};
    LPupilDiametermm_M = dataArray{:, 6};
    RDiaXpx_M = dataArray{:, 7};
    RDiaYpx_M = dataArray{:, 8};
    RPupilDiametermm_M = dataArray{:, 9};
    BPORXpx_M = dataArray{:, 10};
    BPORYpx_M = dataArray{:, 11};
    BEventInfo_M = dataArray{:, 12};  

    % Conversione vettore tempi
    Time_M=Time_M-Time_M(1);
    Time_M=Time_M/1000000;
    Time_M=Time_M+Time_M(2);
    
    % Salvataggio dati.
    save([path(i).folder '\' path(i).name(23:25) '_mobile.mat'],...
        'Time_M', 'LPupilDiametermm_M', 'RPupilDiametermm_M', 'BEventInfo_M')
    
    % Pulizia variabili. 
    clearvars -except path delimiter startRow formatSpec i;
    

end