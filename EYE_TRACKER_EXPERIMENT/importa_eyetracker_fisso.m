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

formatSpec = '%f%C%f%f%f%f%f%f%f%C%C%[^\n\r]';

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
    Time_F = dataArray{:, 1};
    Type_F = dataArray{:, 2};
    Trial_F = dataArray{:, 3};
    LDiaXpx_F = dataArray{:, 4};
    LDiaYpx_F = dataArray{:, 5};
    LPupilDiametermm_F = dataArray{:, 6};
    RDiaXpx_F = dataArray{:, 7};
    RDiaYpx_F = dataArray{:, 8};
    RPupilDiametermm_F = dataArray{:, 9};
    LEventInfo_F = dataArray{:, 10};
    REventInfo_F = dataArray{:, 11};     
    
    % Conversione vettore tempi
    Time_F=Time_F-Time_F(1);
    Time_F=Time_F/1000000;
    Time_F=Time_F+Time_F(2);

    % Salvataggio dati.
    save([path(i).folder '\' path(i).name(16:18) '_fisso.mat'],...
        'Time_F', 'LPupilDiametermm_F', 'RPupilDiametermm_F', 'LEventInfo_F', 'REventInfo_F')
    
    % Pulizia variabili. 
    clearvars -except path delimiter startRow formatSpec i;
    

end