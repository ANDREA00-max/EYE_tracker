
% SCRIPT PER LA SINCRONIZZAZIONE DEI VETTORI ESTRATTI DAL MOBILE E DAL FISSO:
% 
%
%
close all;
clear all;
clc;

% -----------------------------------------------------------------------
% Inizializzazione variabli.
% -----------------------------------------------------------------------

path=dir([pwd,'/Dati/S*']);

for i=1:length(path)

    %CARICO I DATI DEL FISSO 
    load([path(i).folder '/' path(i).name '/' path(i).name '_fisso_ripulito.mat'])
    
    Time_F_Finale=Time_F;
    LPupilDiametermm_F_Finale=LPupilDiametermm_F;
    RPupilDiametermm_F_Finale=RPupilDiametermm_F;
    A=LPupilDiametermm_F_Finale(5000:end);
    B=RPupilDiametermm_F_Finale(5000:end);
    if(i==21 | i ==19 | i==15 | i==14 | i==7 )
        A=A(5000:end);
        B=B(5000:end);
    end
    
    if i==19
        A=A(1000:end);
        B=B(1000:end);
    end
  
    %CARICO I DATI DEL MOBILE 
    load([path(i).folder '/' path(i).name '/' path(i).name '_mobile_aggiornato_ripulito.mat'])
    
    RPupilDiametermm_M_Finale=RPupilDiametermm_M;
    LPupilDiametermm_M_Finale=LPupilDiametermm_M;
    C=LPupilDiametermm_M_Finale(5000:end);
    D=RPupilDiametermm_M_Finale(5000:end);
    if(i==21 | i ==19 | i==15 | i==14 | i==7)
        C=C(5000:end);
        D=D(5000:end);
    end
    
    if i==19
        C=C(1000:end);
        D=D(1000:end);
    end
    
    
    LagSx=finddelay(A, C);
    LagDx=finddelay(B, D);
    
    Lag=min(LagSx,LagDx);
    
    if(i==21 | i ==19 | i==15 | i==7)
        Lag=max(LagSx,LagDx);
    end
    
    LPupilDiametermm_F_Finale=LPupilDiametermm_F_Finale(-Lag:end);
    RPupilDiametermm_F_Finale=RPupilDiametermm_F_Finale(-Lag:end);
    Time_F_Finale=((1:length(LPupilDiametermm_F_Finale))/30)';
    Time_M_Finale=((1:length(LPupilDiametermm_M_Finale))/30)';
    
    
     save([path(i).folder '/' path(i).name '/' path(i).name '_sincronizzati'],...
         'LPupilDiametermm_F_Finale','LPupilDiametermm_M_Finale','RPupilDiametermm_M_Finale','RPupilDiametermm_F_Finale','Time_M_Finale','Time_F_Finale')
   figure;
     plot(LPupilDiametermm_F_Finale)
     hold on
     plot(LPupilDiametermm_M_Finale)
end
