%Script per il calcolo delle matrici di grandi dimensioni che utilizzeremo
%per i test statistici. La normalizzazione dei valori è "commentata"
%perché non sono sicuro si debba eseguire 

clear all;
close all;
clc;

MEDIE_LF=zeros(9,30);
MEDIE_LM=zeros(9,30);
MEDIE_RM=zeros(9,30);
MEDIE_RF=zeros(9,30);

STD_LF=zeros(9,30);
STD_LM=zeros(9,30);
STD_RM=zeros(9,30);
STD_RF=zeros(9,30);


path=dir([pwd,'/Dati/S*']);

for i=1:30
    
    load([path(i).folder '/' path(i).name '/' path(i).name '_medie_max_min.mat'])
    
%     Media_Sx_Fisso=(Media_Sx_Fisso-mean(Media_Sx_Fisso))/std(Media_Sx_Fisso);
%     Media_Sx_Mobile=(Media_Sx_Mobile-mean(Media_Sx_Mobile))/std(Media_Sx_Mobile);
%     Media_Dx_Mobile=(Media_Dx_Mobile-mean(Media_Dx_Mobile))/std(Media_Dx_Mobile);
%     Media_Dx_Fisso=(Media_Dx_Fisso-mean(Media_Dx_Fisso))/std(Media_Dx_Fisso);
    
%     Std_Sx_Fisso=(Std_Sx_Fisso-mean(Std_Sx_Fisso))/std(Std_Sx_Fisso);
%     Std_Sx_Mobile=(Std_Sx_Mobile-mean(Std_Sx_Mobile))/std(Std_Sx_Mobile);
%     Std_Dx_Mobile=(Std_Dx_Mobile-mean(Std_Dx_Mobile))/std(Std_Dx_Mobile);
%     Std_Dx_Fisso=(Std_Dx_Fisso-mean(Std_Dx_Fisso))/std(Std_Dx_Fisso);
    
    MEDIE_LF(:,i)=Media_Sx_Fisso;
    MEDIE_LM(:,i)=Media_Sx_Mobile;
    MEDIE_RM(:,i)=Media_Dx_Mobile;
    MEDIE_RF(:,i)=Media_Dx_Fisso;
    
    STD_LF(:,i)=Std_Sx_Fisso;
    STD_LM(:,i)=Std_Sx_Mobile;
    STD_RM(:,i)=Std_Dx_Mobile;
    STD_RF(:,i)=Std_Dx_Fisso;
    
end


    
 
