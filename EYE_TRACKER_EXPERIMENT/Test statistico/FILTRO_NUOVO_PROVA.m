% ----------------------------------------
% Script per calcolare le potenze spettrali dei segnali divisi nelle 9 zone
% innanzitutto filtraggio dei segnali
% ---------------------------------------------------

clear all;
close all;
clc;

%inizializzo delle variabili

%----VARIABILI FILTRO-----
PBfreq = 1; % Fino a dove il segnale deve essere lasciato passare "indenne"
SBfreq = 1.5; % Dove inizia arresto della banda (da PBfreq a SBfreq abbiamo banda di transizione)
PBripple = 0.5; % Default = 1
SBattenuation = 60; % Default = 60 (-dB)
%-------------------------

%----variabili posizioni arousal costante-----

istante=zeros(9,2);
intervallo_gruppo_Pupilla=istante;
intervallo_gruppo_Tacogramma=istante; 

%sono le matrici che, per ogni soggetto, utilizzerò come appoggio per il
%calcolo delle posizioni dei vettori in cui i valori di arousal e valence
%sono noti per poter eseguire il calcolo della potenza spettrale. Sono
%delle matrici 9x2 perché le zone sono 9, quindi 9 righe, e ci sarà indice
%di inizio e fine di ogni zona (quindi due elementi per riga).

%----------------------------------------------



%1)----FILTRO PASSA BASSO---- ------------------------------------
DesignMethod = 'kaiserwin'; % kaiserwin e' uno dei migliori per filtri FIR (si puo provare anche "equiripple", piï¿½ tradizionale)
fs = 30; % Specificando la freq di campionamento possiamo evitare di ragionare in unita' normalizzate ...
lpFilt = designfilt('lowpassfir','PassbandFrequency',PBfreq, ...
'StopbandFrequency',SBfreq,'PassbandRipple',PBripple, ...
'StopbandAttenuation',SBattenuation,'DesignMethod',DesignMethod,...
'SampleRate',fs);

%1.1) filtro passabasso per tacogramma, a 4 Hz
lpFilt_Taco = designfilt('lowpassfir','PassbandFrequency',PBfreq, ...
'StopbandFrequency',SBfreq,'PassbandRipple',PBripple, ...
'StopbandAttenuation',SBattenuation,'DesignMethod',DesignMethod,...
'SampleRate',4);
%---------------------------------------------------


%2)------------Filtraggio passa altO------------------

PBfreq = 0.025;
order = 3;
hpFilt = designfilt('highpassiir','FilterOrder',order, ...
         'HalfPowerFrequency',PBfreq,'SampleRate',fs);
     
%2.1) Taco
     PBfreq = 0.025;
order = 3;
hpFilt_Taco = designfilt('highpassiir','FilterOrder',order, ...
         'HalfPowerFrequency',PBfreq,'SampleRate',4);
%-------------------------------------------------------

%scelgo la cartelle di riferimento ed inizializzo le variabili

path_eventi=dir([pwd,'/Eventi/esp*']);
path=dir([pwd,'/Dati/S*']);

for i=1:length(path)
    
    load([path(i).folder '/' path(i).name '/' path(i).name '_arousal'])

    sig_PD_LM_prefilt = LPupilDiametermm_M;
    sig_PD_LF_prefilt = LPupilDiametermm_F;
    sig_PD_RM_prefilt = RPupilDiametermm_M;
    sig_PD_RF_prefilt = RPupilDiametermm_F;
    
    %ricampiono a 4 per poter applicare il filtro al segnale
    
    [sig_TACO_prefilt,A]=resample(Taco,tempo_Taco,4);
    [sig_RESP_prefilt,B]=resample(Resp,tempo_Taco,4);
    tempo_Taco=A;
    
%----3)applico i filtri sia passaalto che passabasso-------
    
    sig_LM_lp = filtfilt(lpFilt,sig_PD_LM_prefilt);
    sig_LF_lp = filtfilt(lpFilt,sig_PD_LF_prefilt);
    sig_RM_lp = filtfilt(lpFilt,sig_PD_RM_prefilt);
    sig_RF_lp = filtfilt(lpFilt,sig_PD_RF_prefilt);
    
    sig_RESP_lp=filtfilt(lpFilt_Taco,sig_RESP_prefilt);
    sig_TACO_lp=filtfilt(lpFilt_Taco,sig_TACO_prefilt);
    
    sig_LM_hp = filtfilt(hpFilt,sig_LM_lp);
    sig_LF_hp = filtfilt(hpFilt,sig_LF_lp);
    sig_RM_hp = filtfilt(hpFilt,sig_RM_lp);
    sig_RF_hp = filtfilt(hpFilt,sig_RF_lp);
    
    sig_RESP_hp=filtfilt(hpFilt_Taco,sig_RESP_lp);
    sig_TACO_hp=filtfilt(hpFilt_Taco,sig_TACO_lp);
    
    
%-------------------------------------------------------
    
% 4)ricampionamento a 4 Hz-------------------------------------------------
    
    [y_LF,ty_LF] = resample(sig_LF_hp,tempo_Fisso,4,'pchip');
    [y_RF,ty_RF] = resample(sig_RF_hp,tempo_Fisso,4,'pchip');
    [y_LM,ty_LM] = resample(sig_LM_hp,tempo_Mobile,4,'pchip');
    [y_RM,ty_RM] = resample(sig_RM_hp,tempo_Mobile,4,'pchip');
    [y_TACO,ty_TACO] = resample(sig_TACO_hp,tempo_Taco,4,'pchip');
    [y_RESP,ty_RESP] = resample(sig_RESP_hp,tempo_Taco,4,'pchip'); %i tempi di resp e taco sono uguali!!
%--------------------------------------------------------------------------
    

%5)------per delle convenzioni utilizzo i millisecondi-----
%---------- per comodità eseguo la trasposizione dei tempi-----
     ty_LF=ty_LF'*1000;
    ty_RF=ty_RF'*1000;
    ty_LM=ty_LM'*1000;
    ty_RM=ty_RM'*1000;
    ty_TACO=ty_TACO'*1000;
    ty_RESP=ty_RESP'*1000;  
%---------------------------------------------------

    
%--------------Calcolo degli indici di interesse per analisi della potenza spettrale------------
    
    load([path_eventi(i).folder '/' path_eventi(i).name '/' path_eventi(i).name '_i.mat'])
    
    for k=1:9
        
        j=(k-1)*10;
    
        trova = find(eventi(:,1)>=1+j & eventi(:,1)<=10+j);
        istante(k,1)=eventi(trova(1),2);
        istante(k,2)=eventi(trova(end)+1,2); %mi fermo dove inizia la baseline (+1)
        
        intervallo_gruppo_Pupilla(k,1) = find(ty_LF>istante(k,1)*1000,1,'first');
        intervallo_gruppo_Pupilla(k,2) = find(ty_LF<istante(k,2)*1000,1,'last');
        
    end

     M= intervallo_gruppo_Pupilla;
 
 
 %%tacogramma e respirogramma
 
    for k=1:9
    
        j=(k-1)*10;
    
        trova = find(eventi(:,1)>=1+j & eventi(:,1)<=10+j);
        istante(k,1)=eventi(trova(1),2);
        istante(k,2)=eventi(trova(end)+1,2); %mi fermo dove inizia la baseline (+1)
        
    
        intervallo_gruppo_Tacogramma(k,1) = find(ty_TACO>istante(k,1)*1000,1,'first');
        intervallo_gruppo_Tacogramma(k,2) = find(ty_TACO<istante(k,2)*1000,1,'last');
    
    end
    T=intervallo_gruppo_Tacogramma;
    
%----------------------------------------------------------------------------------
    

    SEGNALE_LF = [ty_LF,y_LF];
    SEGNALE_RF = [ty_RF,y_RF];
    SEGNALE_LM = [ty_LM,y_LM];
    SEGNALE_RM = [ty_RM,y_RM];
    SEGNALE_TACO = [ty_TACO,y_TACO];
    SEGNALE_RESP = [ty_RESP, y_RESP];
   
    %salvo gli outuput della powerunivariateAR
    
     [BandPower1,stats1] = powerUnivariateAR(SEGNALE_LF(M(1,1):M(1,2),:),4,16);
     [BandPower2,stats2] = powerUnivariateAR(SEGNALE_LF(M(2,1):M(2,2),:),4,16);
     [BandPower3,stats3] = powerUnivariateAR(SEGNALE_LF(M(3,1):M(3,2),:),4,16);
     [BandPower4,stats4] = powerUnivariateAR(SEGNALE_LF(M(4,1):M(4,2),:),4,16);
     [BandPower5,stats5] = powerUnivariateAR(SEGNALE_LF(M(5,1):M(5,2),:),4,16);
     [BandPower6,stats6] = powerUnivariateAR(SEGNALE_LF(M(6,1):M(6,2),:),4,16);
     [BandPower7,stats7] = powerUnivariateAR(SEGNALE_LF(M(7,1):M(7,2),:),4,16);
     [BandPower8,stats8] = powerUnivariateAR(SEGNALE_LF(M(8,1):M(8,2),:),4,16);
     [BandPower9,stats9] = powerUnivariateAR(SEGNALE_LF(M(9,1):M(9,2),:),4,16);  
     
     BandPowerLF=[BandPower1'; BandPower2'; BandPower3'; BandPower4'; BandPower5'; BandPower6'; BandPower7'; BandPower8'; BandPower9'];
     statsLF=[stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9];
     
    [BandPower1,stats1] = powerUnivariateAR(SEGNALE_RF(M(1,1):M(1,2),:),4,16);
    [BandPower2,stats2] = powerUnivariateAR(SEGNALE_RF(M(2,1):M(2,2),:),4,16);
    [BandPower3,stats3] = powerUnivariateAR(SEGNALE_RF(M(3,1):M(3,2),:),4,16);
    [BandPower4,stats4] = powerUnivariateAR(SEGNALE_RF(M(4,1):M(4,2),:),4,16);
    [BandPower5,stats5] = powerUnivariateAR(SEGNALE_RF(M(5,1):M(5,2),:),4,16);
    [BandPower6,stats6] = powerUnivariateAR(SEGNALE_RF(M(6,1):M(6,2),:),4,16);
    [BandPower7,stats7] = powerUnivariateAR(SEGNALE_RF(M(7,1):M(7,2),:),4,16);
    [BandPower8,stats8] = powerUnivariateAR(SEGNALE_RF(M(8,1):M(8,2),:),4,16);
    [BandPower9,stats9] = powerUnivariateAR(SEGNALE_RF(M(9,1):M(9,2),:),4,16);
    
    BandPowerRF=[BandPower1'; BandPower2'; BandPower3'; BandPower4'; BandPower5'; BandPower6'; BandPower7'; BandPower8'; BandPower9'];
    statsRF=[stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9];
     
     
    [BandPower1,stats1] = powerUnivariateAR(SEGNALE_LM(M(1,1):M(1,2),:),4,16);
    [BandPower2,stats2] = powerUnivariateAR(SEGNALE_LM(M(2,1):M(2,2),:),4,16);
    [BandPower3,stats3] = powerUnivariateAR(SEGNALE_LM(M(3,1):M(3,2),:),4,16);
    [BandPower4,stats4] = powerUnivariateAR(SEGNALE_LM(M(4,1):M(4,2),:),4,16);
    [BandPower5,stats5] = powerUnivariateAR(SEGNALE_LM(M(5,1):M(5,2),:),4,16);
    [BandPower6,stats6] = powerUnivariateAR(SEGNALE_LM(M(6,1):M(6,2),:),4,16);
    [BandPower7,stats7] = powerUnivariateAR(SEGNALE_LM(M(7,1):M(7,2),:),4,16);
    [BandPower8,stats8] = powerUnivariateAR(SEGNALE_LM(M(8,1):M(8,2),:),4,16);
    [BandPower9,stats9] = powerUnivariateAR(SEGNALE_LM(M(9,1):M(9,2),:),4,16);
    
     BandPowerLM=[BandPower1'; BandPower2'; BandPower3'; BandPower4'; BandPower5'; BandPower6'; BandPower7'; BandPower8'; BandPower9'];
    statsLM=[stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9];
    
    [BandPower1,stats1] = powerUnivariateAR(SEGNALE_RM(M(1,1):M(1,2),:),4,16);
    [BandPower2,stats2] = powerUnivariateAR(SEGNALE_RM(M(2,1):M(2,2),:),4,16);
    [BandPower3,stats3] = powerUnivariateAR(SEGNALE_RM(M(3,1):M(3,2),:),4,16);
    [BandPower4,stats4] = powerUnivariateAR(SEGNALE_RM(M(4,1):M(4,2),:),4,16);
    [BandPower5,stats5] = powerUnivariateAR(SEGNALE_RM(M(5,1):M(5,2),:),4,16);
    [BandPower6,stats6] = powerUnivariateAR(SEGNALE_RM(M(6,1):M(6,2),:),4,16);
    [BandPower7,stats7] = powerUnivariateAR(SEGNALE_RM(M(7,1):M(7,2),:),4,16);
    [BandPower8,stats8] = powerUnivariateAR(SEGNALE_RM(M(8,1):M(8,2),:),4,16);
    [BandPower9,stats9] = powerUnivariateAR(SEGNALE_RM(M(9,1):M(9,2),:),4,16);
     
    BandPowerRM=[BandPower1'; BandPower2'; BandPower3'; BandPower4'; BandPower5'; BandPower6'; BandPower7'; BandPower8'; BandPower9'];
   statsRM=[stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9];
    
    [BandPower1,stats1] = powerUnivariateAR(SEGNALE_TACO(T(1,1):T(1,2),:),4,16);
    [BandPower2,stats2] = powerUnivariateAR(SEGNALE_TACO(T(2,1):T(2,2),:),4,16);
    [BandPower3,stats3] = powerUnivariateAR(SEGNALE_TACO(T(3,1):T(3,2),:),4,16);
    [BandPower4,stats4] = powerUnivariateAR(SEGNALE_TACO(T(4,1):T(4,2),:),4,16);
    [BandPower5,stats5] = powerUnivariateAR(SEGNALE_TACO(T(5,1):T(5,2),:),4,16);
    [BandPower6,stats6] = powerUnivariateAR(SEGNALE_TACO(T(6,1):T(6,2),:),4,16);
    [BandPower7,stats7] = powerUnivariateAR(SEGNALE_TACO(T(7,1):T(7,2),:),4,16);
    [BandPower8,stats8] = powerUnivariateAR(SEGNALE_TACO(T(8,1):T(8,2),:),4,16);
    [BandPower9,stats9] = powerUnivariateAR(SEGNALE_TACO(T(9,1):T(9,2),:),4,16);
    
    BandPowerTACO=[BandPower1'; BandPower2'; BandPower3'; BandPower4'; BandPower5'; BandPower6'; BandPower7'; BandPower8'; BandPower9'];
    statsTACO=[stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9];
    
    [BandPower1,stats1] = powerUnivariateAR(SEGNALE_RESP(T(1,1):T(1,2),:),4,16);
    [BandPower2,stats2] = powerUnivariateAR(SEGNALE_RESP(T(2,1):T(2,2),:),4,16);
    [BandPower3,stats3] = powerUnivariateAR(SEGNALE_RESP(T(3,1):T(3,2),:),4,16);
    [BandPower4,stats4] = powerUnivariateAR(SEGNALE_RESP(T(4,1):T(4,2),:),4,16);
    [BandPower5,stats5] = powerUnivariateAR(SEGNALE_RESP(T(5,1):T(5,2),:),4,16);
    [BandPower6,stats6] = powerUnivariateAR(SEGNALE_RESP(T(6,1):T(6,2),:),4,16);
    [BandPower7,stats7] = powerUnivariateAR(SEGNALE_RESP(T(7,1):T(7,2),:),4,16);
    [BandPower8,stats8] = powerUnivariateAR(SEGNALE_RESP(T(8,1):T(8,2),:),4,16);
    [BandPower9,stats9] = powerUnivariateAR(SEGNALE_RESP(T(9,1):T(9,2),:),4,16);
   
    BandPowerRESP=[BandPower1'; BandPower2'; BandPower3'; BandPower4'; BandPower5'; BandPower6'; BandPower7'; BandPower8'; BandPower9'];
    statsRESP=[stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9];


%------------------------salvo tutti i risultati---------------------------------

save([path(i).folder '/' path(i).name '/' path(i).name '_arousal_con_psd_PROVA'],...
         'LPupilDiametermm_F','RPupilDiametermm_F','RPupilDiametermm_M','LPupilDiametermm_M',...
         'Taco','Resp',...
         'tempo_Mobile','tempo_Fisso','tempo_Taco',...
         'BandPowerRESP','BandPowerTACO','BandPowerLM','BandPowerRM','BandPowerRF','BandPowerLF',...
         'statsRESP','statsTACO','statsLM','statsRM','statsRF','statsLF','M','T')
     
   
end