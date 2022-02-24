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

     elimina_occhi_chiusi = find(ty_LF>=eventi(3,2)*1000,1,'first');
     M= intervallo_gruppo_Pupilla - elimina_occhi_chiusi;
 
 
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
   
    %da qui in poi grafico i risultati
    
for j=1:9
    
     [BandPower,stats] = powerUnivariateAR(SEGNALE_LF(M(j,1):M(j,2),:),4,6:30);
    % figure;
   freq(:,j)=stats.PSDfreqs;
     psd_LF(:,j)=stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1});
%     plot(freq_LP(:,j),psd_LP(:,j))
    % title('pupilla sinistra')
  
end
% %      
    for j=1:9
    
    [BandPower,stats] = powerUnivariateAR(SEGNALE_RF(M(j,1):M(j,2),:),4,6:30);
    psd_RF(:,j)=stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1});
    
%     figure;
%     plot(stats.PSDfreqs,stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1}))
%    title('pupilla destra'); 
    end
%     
    for j=1:9
%     
    [BandPower,stats] = powerUnivariateAR(SEGNALE_LM(M(j,1):M(j,2),:),4,6:30);;
     psd_LM(:,j)=stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1});
%     figure;
%     plot(stats.PSDfreqs,stats.cPSDtotTimeFrames{1})
%     
    end
% 
    for j=1:9
%     
    [BandPower,stats] = powerUnivariateAR(SEGNALE_RM(M(j,1):M(j,2),:),4,6:30);
%     figure;
%     plot(stats.PSDfreqs,stats.cPSDtotTimeFrames{1})
     psd_RM(:,j)=stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1});
    end
%    
%     
%     
% 
    for j=1:9
%     
    [BandPower,stats] = powerUnivariateAR(SEGNALE_TACO(T(j,1):T(j,2),:),4,6:30);
%     figure;
    psd_TACO(:,j)=stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1});
%     plot(freq_T(:,j),psd_T(:,j))
%      title('cuore');
    end
    
    for j=1:9
    
    [BandPower,stats] = powerUnivariateAR(SEGNALE_RESP(T(j,1):T(j,2),:),4,6:30);
%     figure;
     psd_RESP(:,j)=stats.cPSDtotTimeFrames{1}/max(stats.cPSDtotTimeFrames{1});
%     plot(freq_RESP(:,j),psd_RESP(:,j))
%     title('respiro')
     
%------------------------salvo tutti i risultati---------------------------------

 save([path(i).folder '/' path(i).name '/' path(i).name '_arousal'],...
        'LPupilDiametermm_F','RPupilDiametermm_F','RPupilDiametermm_M','LPupilDiametermm_M',...
        'Taco','Resp',...
        'tempo_Mobile','tempo_Fisso','tempo_Taco',...
        'freq','psd_RESP','psd_TACO','psd_LF','psd_RF','psd_RM','psd_LM')
     
  end
   
end