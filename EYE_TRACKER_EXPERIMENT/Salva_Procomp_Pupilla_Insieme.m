
%--------------------------------------------------------------------
% script per salvare i valori procomp insieme a quelli della pupilla
%----------------------------------------------------------------------
close all;
clear;
clc;

% carico il path
% inizializzo delle variabili


path_dati=dir([pwd,'/Dati/S*']);

for i=1:length(path_dati)

    load([path_dati(i).folder '/' path_dati(i).name '/' path_dati(i).name '_ProComp.mat'])
    load([path_dati(i).folder '/' path_dati(i).name '/' path_dati(i).name '_sincronizzazione_finale.mat'])
    
    tempo_Taco=sig_tacogramma(:,1)';
    Taco=sig_tacogramma(:,2);
    Resp=sig_respirogramma(:,2);
    
      save([path_dati(i).folder '/' path_dati(i).name '/' path_dati(i).name '_arousal'],...
         'LPupilDiametermm_F','RPupilDiametermm_F','RPupilDiametermm_M','LPupilDiametermm_M',...
         'Taco','Resp',...
         'tempo_Mobile','tempo_Fisso','tempo_Taco')
end

