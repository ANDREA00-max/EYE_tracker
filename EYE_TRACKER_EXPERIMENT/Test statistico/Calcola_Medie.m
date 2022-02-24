%-----------------------------------------------------------------------
%% Calcolo medie dei valori di diametro pupillare a seconda dell'arousal
%-----------------------------------------------------------------------
%All'interno di ogni matrice, che sarà una 6x2, avrò tutti i valori
%ordinati come segue:
% sulla colonna sinistra tutti i valori relativi alla pupilla sinistra
% a destra tutti i valori della destra
%i primi tre elementi per colonna relativi al fisso
%gli ultimi tre relativi al mobile
%prima e quarta riga basso arousal
%seconda e quinta riga medio arousal
%terza e sesta riga alto arousal

%Nei vettori Media_Sx_Fisso, Media_Sx_Mobile....
%Std_Sx_Fisso, Std_Sx_Mobile....
%Ci sono i 9 valori relativi alle zone ad arousal costante (Che sono
%appunto 9). Ad ogni soggetto e per ogni misurazione è associato un vettore con 9 valori, nello
%script successivo ("MEDIESTD") unirò questi vettori in matrici più grandi
%(che saranno in tutto 8: 4 per le medie e 4 per le deviazioni standard
%(LF,RF,LM,RM).

clear all;
close all;
clc

path=dir([pwd,'/Dati/S*']);

Media_Sx_Fisso=[];
Media_Sx_Mobile=[];
Media_Dx_Fisso=[];
Media_Dx_Mobile=[];
Std_Sx_Fisso=[];
Std_Sx_Mobile=[];
Std_Dx_Fisso=[];
Std_Dx_Mobile=[];

for k=1:length(path)
    clear Media_Sx_Fisso;
    clear Media_Sx_Mobile;
    clear Media_Dx_Fisso;
    clear Media_Dx_Mobile;
    clear Std_Sx_Fisso;
    clear Std_Sx_Mobile;
    clear Std_Dx_Fisso;
    clear Std_Dx_Mobile;
    
    Media_Sx_Fisso=[];
    Media_Sx_Mobile=[];
    Media_Dx_Fisso=[];
    Media_Dx_Mobile=[];
    Std_Sx_Fisso=[];
    Std_Sx_Mobile=[];
    Std_Dx_Fisso=[];
    Std_Dx_Mobile=[];
    
    load([path(k).folder '/' path(k).name '/' path(k).name '_arousal.mat'])
    for j=1:3
        for i = (j*3-2):(j*3)
            medLF(mod(i,3)+1) = mean (LPupilDiametermm_F((M1(i,1):M1(i,2))));
            pit(mod(i,3)+1) =  max (LPupilDiametermm_F((M1(i,1):M1(i,2))));
            pie(mod(i,3)+1) =  min (LPupilDiametermm_F((M1(i,1):M1(i,2))));
            stdLF(mod(i,3)+1) =  std (LPupilDiametermm_F((M1(i,1):M1(i,2))));
        end
        Media(j,1) = mean (medLF);
        Minimi(j,1) = min (pie);
        Massimi (j,1) = max (pit);
        STD (j,1) = mean (stdLF);
        Media_Sx_Fisso=[Media_Sx_Fisso; medLF(1); medLF(2); medLF(3)];
        Std_Sx_Fisso=[Std_Sx_Fisso; stdLF(1); stdLF(2); stdLF(3)];
        
        for i = (j*3-2):(j*3)

            medRF(mod(i,3)+1) = mean (RPupilDiametermm_F((M1(i,1):M1(i,2))));
            pit(mod(i,3)+1) =  max (RPupilDiametermm_F((M1(i,1):M1(i,2))));
            pie(mod(i,3)+1) =  min (RPupilDiametermm_F((M1(i,1):M1(i,2))));
            stdRF(mod(i,3)+1) =  std (RPupilDiametermm_F((M1(i,1):M1(i,2))));
            
        end
        
        Media(j,2) = mean (medRF);
        Minimi(j,2) = min (pie);
        Massimi (j,2) = max (pit);
        STD (j,2) = mean (stdRF);
        Media_Dx_Fisso=[Media_Dx_Fisso; medRF(1); medRF(2); medRF(3)];
        Std_Dx_Fisso=[Std_Dx_Fisso; stdRF(1); stdRF(2); stdRF(3)];
        
        
        for i = (j*3-2):(j*3)

            medLM(mod(i,3)+1) = mean (LPupilDiametermm_M((M1(i,1):M1(i,2))));
            pit(mod(i,3)+1) =  max (LPupilDiametermm_M((M1(i,1):M1(i,2))));
            pie(mod(i,3)+1) =  min (LPupilDiametermm_M((M1(i,1):M1(i,2))));
            stdLM(mod(i,3)+1) =  std (LPupilDiametermm_M((M1(i,1):M1(i,2))));
            
        end
        
        Media(j+3,1) = mean (medLM);
        Minimi(j+3,1) = min (pie);
        Massimi (j+3,1) = max (pit);
        STD (j+3,1) = mean (stdLM);
        Media_Sx_Mobile=[Media_Sx_Mobile; medLM(1); medLM(2); medLM(3)];
        Std_Sx_Mobile=[Std_Sx_Mobile; stdLM(1); stdLM(2); stdLM(3)];
        
        for i = (j*3-2):(j*3)
  
            medRM(mod(i,3)+1) =  mean (RPupilDiametermm_M((M1(i,1):M1(i,2))));
            pit(mod(i,3)+1) =  max (RPupilDiametermm_M((M1(i,1):M1(i,2))));
            pie(mod(i,3)+1) =  min (RPupilDiametermm_M((M1(i,1):M1(i,2))));
            stdRM(mod(i,3)+1) =  std (RPupilDiametermm_M((M1(i,1):M1(i,2))));
        end
        
        Media(j+3,2) = mean (medRM);
        Minimi(j+3,2) = min (pie);
        Massimi (j+3,2) = max (pit);
        STD (j+3,2) = mean (stdRM);
        Media_Dx_Mobile=[Media_Dx_Mobile; medRM(1); medRM(2); medRM(3)];
        Std_Dx_Mobile=[Std_Dx_Mobile; stdRM(1); stdRM(2); stdRM(3)];
        
        
    end
    
     save([path(k).folder '/' path(k).name '/' path(k).name '_medie_max_min.mat'],...
        'Media', 'Massimi', 'Minimi', 'STD', 'Media_Sx_Mobile', 'Media_Dx_Fisso','Media_Dx_Mobile', 'Media_Sx_Fisso',...
        'Std_Dx_Mobile','Std_Sx_Mobile','Std_Dx_Fisso','Std_Sx_Fisso');
end
