%SCRIPT AGGIORNATO PER IL CALCOLO DEGLI INDICI TEMPORALI CORRISPONDENTI ALLE ZONE AD AROUSAL COSTANTE 




%%PER ANALISI STATISTICA, IN TEORIA, I DATI TACO E RESP NON CI SERVONO 
 %ANALISI STATISTICA= MEDIA, MAX, MIN, STD delle pupille nelle zone ad
 %arousal costante.  il codice va in parte modificato per ottenere anche
 %gli indici delle zone di BASELINE e di LOW AROUSAL, MEDIUM AROUSAL, HIGH
 %AROUSAL. 
 %prova a vedere come si fanno i boxplot....

path_eventi=dir([pwd,'/Eventi/esp*']);
path_dati=dir([pwd,'/Dati/S*']);

for i=1:1
    
    load([path_eventi(i).folder '/' path_eventi(i).name '/' path_eventi(i).name '_i.mat'])
        load([path_dati(i).folder '/' path_dati(i).name '/' path_dati(i).name '_arousal.mat'])
    
    for k=1:9
        j=(k-1)*10;
    
        trova = find(eventi(:,1)>=1+j & eventi(:,1)<=10+j);
        istante(k,1)=eventi(trova(1),2);
        istante(k,2)=eventi(trova(end)+1,2); %mi fermo dove inizia la baseline (+1)
        
        intervallo_gruppo_Pupilla(k,1) = find(tempo_Fisso>istante(k,1),1,'first');
        intervallo_gruppo_Pupilla(k,2) = find(tempo_Fisso<istante(k,2),1,'last');
        
    end

     elimina_occhi_chiusi = find(tempo_Fisso>=eventi(3,2)*1000,1,'first');
     M= intervallo_gruppo_Pupilla - elimina_occhi_chiusi;
 
 
 %%tacogramma e respirogramma
 
    for k=1:9
    
        j=(k-1)*10;
    
        trova = find(eventi(:,1)>=1+j & eventi(:,1)<=10+j);
        istante(k,1)=eventi(trova(1),2);
        istante(k,2)=eventi(trova(end)+1,2); %mi fermo dove inizia la baseline (+1)
        
    
        intervallo_gruppo_Tacogramma(k,1) = find(tempo_Taco>istante(k,1),1,'first');
        intervallo_gruppo_Tacogramma(k,2) = find(tempo_Taco<istante(k,2),1,'last');
    
    end
    T=intervallo_gruppo_Tacogramma;
    
end

