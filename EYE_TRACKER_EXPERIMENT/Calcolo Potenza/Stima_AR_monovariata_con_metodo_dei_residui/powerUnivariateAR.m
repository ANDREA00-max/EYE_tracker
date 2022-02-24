function [BandPower,stats] = powerUnivariateAR(taco,fs,order,freqbands,epochLen,overlapPerc,vEventi,residWhiteTestOpt)

%%%% DESCRIZIONE:
% Calcolo della potenza con stima AR monovariata e metodo dei residui (o
% anche metodo dell'area, decommentando la parte di codice dedicata).
% Riferimenti: Task Force 1996
%              Kamath 2013
%
%%%% USO:
% [BandPower,stats] = powerUnivariateAR(taco,fs,order,freqbands,epochLen,overlapPerc,vEventi,residWhiteTestOpt)
%
%%%% INPUTS:
% Solo i PRIMI DUE SONO OBBLIGATORI, gli altri sono tutti OPZIONALIsig_tacogramma (si possono
% non dichiarare o dichiarare come vettori vuoti per usare i valori di
% default).
%
% taco = matrice del segnale HRV (serie temporale degli intervalli RR),
%        avente nella prima colonna gli istanti di occorrenza del picco R
%        in MILLISECONDI e nella seconda colonna le distanze inter-battito
%        in MILLISECONDI (da un punto di vista di pre-processing e
%        processing del segnale non cambierebbe assolutamente nulla: e'
%        solo questione di convenzione, dalla quale dipende il calcolo del
%        valore di RSA finale [rispettiamo la convenzione di Porges, che
%        utilizza i MILLISECONDI]).
%        Il tacogramma in ingresso DEVE essere GIA' stato RICAMPIONATO e
%        FILTRATO ("detrendato") nel modo desiderato.
%
% fs = frequenza in Hz a cui e' stato ricampionato il taco in input.
%
% order = ordine del modello AR. Puo' essere:
%         vettore     --> Riga o colonna. Stima automatica dell'ordine
%                         "ottimo" per il taco dato in ingresso (ordine piu'
%                         frequente tra quelli che minimizzano AIC su tutti
%                         i timeFrames considerati, tra tutte le fasi della
%                         prova), tra gli ordini passati in ingresso con
%                         questo vettore.
%         num. intero --> la stima del modello AR verra' effettuata con
%                         l'ordine richiesto.
%
% freqbands = matrice contenente gli estremi iniziale e finale (in Hz)
%             della/e banda/e di frequenza in cui calcolare la potenza del
%             taco: una riga per ogni banda di interesse e due colonne per
%             specificare gli estremi minimo (incluso) e massimo (escluso)
%             di ogni banda.
%             Per la stima della RSA, dipende dalle frequenze in cui ci si
%             aspetta di trovare effetto del respiro sul taco nel gruppo di
%             soggetti analizzati. Esempi:
%             Standard per gli adulti --> [0.12-0.40] (Lewis, 2012);
%             "Young children" --> [0.24,1.04] (Miller, 2013).
%
% epochLen = durata degli spezzoni di segnale (in MILLISECONDI, per coerenza
%            con l'unita' di misura del taco) su cui calcolare il valore di
%            potenza in ogni banda richiesta. I valori ottenuti per ogni
%            epoca saranno mediati tra loro per ottenere il valore di
%            potenza finale relativo a ciascuna fase.
%
% overlapPerc = percentuale di overlap tra epoche successive di segnale
%               considerate per il calcolo della potenza. Puo' variare tra
%               0 (nessun overlap) e 99 (overlap quasi completo).
%
% vEventi = vettore contenente l'istante temporale (in MILLISECONDI, per
%           coerenza con l'unita' di misura del taco) di inizio e fine di
%           ciascuna fase della prova da analizzare. Con la seguente
%           struttura:
%           1. INIZIO Fase1
%           2. FINE Fase1
%           3. INIZIO Fase2
%           4. FINE Fase2
%           5. ...
%           Se non passato in ingresso o se passato come vettore vuoto,
%           viene assunto un unico evento all'interno dello spezzone di
%           segnale HRV passato alla funzione.
%
% residWhiteTestOpt = struct contenente i parametri per il funzionamento
%                     del test sui residui. Campi ammessi nella struct:
%                     method --> specifica il metodo con cui deve essere
%                                effettuato test di bianchezza dei residui
%                                (vedi funzione "andersonWhiteTest.m" per
%                                le opzioni disponibili).
%                     alpha  --> coefficiente di significativita' con cui
%                                deve essere svolto test di bianchezza.
%
%%%% OUTPUTS:
% BandPower = valore di potenza calcolato, nelle bande di interesse, con il
%             metodo di stima AR monovariata per il tacogramma fornito in
%             ingresso. Avra' unita' di misura: [ms^2] Per il calcolo della
%             RSA, per favorire confronto con metodo di Porges, si
%             consiglia di applicare logaritmo all'esterno di questa
%             funzione: log(BandPower) --> [ln(ms^2)].
%             NB: Un confronto DIRETTO con RSA ottenuta con metodo di Porges non
%                 sara' comunque possibile, neanche a parit� di "epochsLen", dato
%                 il diverso calcolo della potenza del taco effettuato dai due
%                 metodi.
%
% stats = struct contenente cell-arrays con parametri diagnostici vari.
%         Ogni cella del cell-array si riferisce ad ogni fase della prova
%         (definita da ogni coppia consecutiva di eventi "vEventi"), mentre
%         ogni elemento contenuto nei vettori presenti in ogni cella si
%         riferisce ad una singola epoca di segnale (di lunghezza epochLen)
%         analizzata per la singola fase.
%         Precisamente, la struct contiene i seguenti campi (cell array):
%         cNonStatTimeFrames --> cell array contenente una matrice che, per
%                                ogni riga, specifica inizio e fine (in ms)
%                                delle epoche di segnale scartate per
%                                violazione dell'ipotesi di stazionarieta'.
%         cStatTimeFrames    --> come sopra, ma per ogni riga contiene
%                                istante di inizio e fine delle epoche di
%                                segnale stazionarie (quelle effettivamente
%                                utilizzate per calcolare RSA)
%         cMeanTacoTimeFrames--> cell array contenente, per ogni fase
%                                (cella), un vettore con il valore medio di
%                                distanza RR rilevata sul tacogramma in
%                                input.
%         cStdTacoTimeFrames --> come sopra, ma riferito alla deviazione
%                                standard.
%         cCImeanTimeFrames  --> cell array contenente, per ogni fase
%                                (cella), un vettore di 2 righe con
%                                limite superiore (riga 1) e inferiore
%                                (riga 2) dell'intervgallo di confidenza
%                                sulla media calcolato per la singola fase.
%         cCIstdTimeFrames   --> come sopra, ma riferito alla deviazione
%                                standard.
%         cOrdOptTimeFrames --> cell array contenente, per ogni fase, un
%                               vettore con l'ordine ottimo individuato con
%                               il criterio di Aikake per ciascuna epoca di
%                               segnale analizzata.
%         cPosOrdOptTimeFrames--> cell array contenente, per ogni fase, un
%                               vettore con la minima cifra AIC scelta (la
%                               prima, o la seconda, o la terza, ecc.),
%                               combinando il criterio di Aikake con la
%                               verifica dell'ipotesi di bianchezza dei
%                               residui. Una cella per ciascuna epoca di
%                               segnale analizzata.
%         cBandPowerTimeFrames--> cell array contenente, per ogni fase, i
%                                 valori di potenza stimati per ogni epoca
%                                 (colonna della matrice) e per ogni banda
%                                 richiesta (riga).
%         cPowerPolesTimeFrames-->cell array contenente, per ogni fase, la
%                                 potenza stimata per ogni polo (riga della
%                                 matrice) e per ogni epoca (colonna).
%         cFreqPolesTimeFrames--> cell array contenente, per ogni fase, la
%                                 frequenza associata a ciascun polo
%                                 contenuto nel cell array precedente.
%         cPSDtotTimeFrames --> cell array contenente, per ogni fase, la
%                               PSD complessiva stimata sommando i
%                               contributi di tutti i poli.
%         PSDfreqs --> asse delle frequenze da usare per visualizzare la PSDtot
%                      relativa a ciascuna fase della prova.
%         cAtestTimeFrames_h--> cell array contenente, per ogni fase, un
%                               vettore con indicazione sul superamento del
%                               test di bianchezza per ogni epoca di
%                               segnale analizzata (si riferisce al modello
%                               AR stimato con il singolo ordine ottimo
%                               scelto per analizzare l'intera traccia).
%                               Il test di bianchezza e' effettuato
%                               considerando il singolo ordine ottimo del
%                               modello scelto per l'intera traccia!
%         cAtestTimeFrames_p--> cell array contenente, per ogni fase, un
%                               vettore con indicazione del p-value
%                               calcolato per il test di bianchezza per
%                               ogni epoca di segnale analizzata (vedi
%                               parentesi per campo struct precedente).
%                               Il test di bianchezza e' effettuato
%                               considerando il singolo ordine ottimo del
%                               modello scelto per l'intera traccia!
%         finalOrder        --> scalare contente l'ordine utilizzato per
%                               stimare la RSA (se input "order" e'
%                               scalare, il valore e' lo stesso dato in
%                               input, altrimenti sara' il valore ottimo
%                               scelto automaticamente tra i valori del
%                               vettore "order" dato in input).
%
%%%%% INFORMAZIONI:
% Autore: Pierluigi Reali, PhD candidate
%         Dipartimento di Elettronica Informazione e Bioingegneria (DEIB)
%         Politecnico di Milano
% Anno: 2018
% Email: pierluigi.reali@polimi.it
%
%**************************************************************************

%0. Controlli e inizializzazioni parametri di default
orderDef = 1:20; % Stima automatica dell'ordine ottimo per il singolo taco
                 % passato in ingresso facendolo variare tra 1 e 20
                 % Riferimenti: Suggerimento Giulia a tesista Adriana (2018): 1-15. 
                 %              Task Force 1996: 8-20.
bandsDef = [0.04,0.15; 0.15,0.4]; % Valida per soggetti adulti.
                                 % Riferimento: Task Force, 1996
epochLenDef = 2 * 60 * 10^3; % Riferimento: Task Force, 1996:
                             % Sono necessari almeno 2 minuti per indagare
                             % le LF. Per stimare la RSA, pero', non sarebbe
                             % necessario indagare le LF ma soltanto la
                             % banda HF, quindi la durata di ogni epoca
                             % considerata si potrebbe tranquillamente
                             % restringere.
overlapPercDef = 50; % Riferimento: Verificare la percentuale consigliata su
                     % qualche fonte (es. Kamath 2013)

if nargin <3 || isempty(order)
    order = orderDef;
elseif ~isnumeric(order) || ~isvector(order) || ... %"isvector(X)"=1 anche se X e' scalare
       ~isequal(order,round(order)) %Quest'ultimo verifica che "order" sia un array di interi
                                    %(precisamente: single/double senza parte decimale)
    error('powerUnivariateAR:orderTypeChk','Parametro order non valido');
end
if nargin <4 || isempty(freqbands)
    freqbands = bandsDef;
elseif ~isnumeric(freqbands) || size(freqbands,2)~=2
    error('powerUnivariateAR:bandChk','Parametro freqbands non valido');
end
if nargin <5 || isempty(epochLen)
    epochLen = epochLenDef;
end
if nargin <6 || isempty(overlapPerc)
    overlapPerc = overlapPercDef;
end
if nargin <7
    vEventi = [];
elseif ~(isvector(vEventi) && isnumeric(vEventi))
    error('powerUnivariateAR:vEventiTypeChk','Parametro vEventi non valido');
elseif rem(length(vEventi),2)~=0 %Inoltre, verifichiamo che la sua lunghezza sia pari
    error('powerUnivariateAR:vEventiTypeChk','La dimensione del parametro vEventi deve essere pari');
end
if nargin <8 || isempty(residWhiteTestOpt)
    residWhiteTestOpt.method = 'adtest';
    residWhiteTestOpt.alpha = 0.05;
end

%--------------------------------------------------------------------------
% 1. Calcoliamo i time frames degli spezzoni di segnale da considerare:
%    viene restituito un cell array "timeFrames" contenente, in ciascuna
%    cella, una matrice con gli istanti iniziale e finale (in MILLISECONDI)
%    di tutte le epoche da considerare per una singola fase della prova
%    (coppia di eventi consevutivi nel vettore "vEventi"). Tali time frames
%    verranno considerati sia per trovare l'ordine ottimo del modello AR
%    stimato dal taco dato in input (se order='auto'), sia per effettuare
%    la stima finale della RSA per ogni fase della prova.

% Overlap percentuale delle epoche, rispetto alla durata di ciascuna di
% esse (epochLen):
overlap = epochLen*(overlapPerc/100);

% Se il vettore degli eventi e' stato passato in input vuoto (oppure se non
% e' stato proprio passato), creo un vettore fittizio che ha come primo
% evento (inizio fase) l'istante iniziale del tacogramma e come secondo evento
% (termine fase) l'istante finale dello stesso:
if isempty(vEventi)
    vEventi=[taco(1,1);taco(end,1)];
end

% Calcoliamo i time frames degli spezzoni di segnale da considerare,
% definiti per ciascuna fase della prova:
nfasi = length(vEventi)/2;
timeFrames = cell([1,nfasi]);
c=1;
for k=1:2:length(vEventi)-1 %Ciclo sugli eventi del vettore, di due in due
    tStart = vEventi(k);
    tEnd = tStart+epochLen;
    timeFrames{c} = [];
    while tEnd <= vEventi(k+1)
        timeFrames{c} = [timeFrames{c};tStart,tEnd];
        tStart = tEnd-overlap;
        tEnd = tStart+epochLen;
    end
    % Se la fase processata non avesse almeno durata pari a epochLen
    % (cioe', fosse piu' corta), allora definisco un'unico timeFrame per
    % quella fase, di durata pari alla sua lunghezza (mando un warning per
    % avvisare):
    if isempty(timeFrames{c})
        warning('powerUnivariateAR:tooShortEvent',...
                'One of the timeFrames has minor length than expected (difference = %f seconds)',...
                (tEnd-vEventi(k+1))/1000);
        timeFrames{c} = [tStart,vEventi(k+1)];
    end
    c=c+1;
end

%--------------------------------------------------------------------------
% % 2. Verifico la stazionariet� del taco nei time frames individuati al passo
% %    precedente: scarto quei time frames per i quali tale stazionarieta'
% %    non risulta verificata con il metodo automatico scelto
% 
% % Parametri per verifica stazionarieta'
% stationarityChkOptions.method = 'bootstrapCircular';
% stationarityChkOptions.alpha = 0.05;
% stationarityChkOptions.bootReps = 50000;
% stationarityChkOptions.blockLen = epochLen/(10^3) * fs;
% stationarityChkOptions.bootciType = 'norm';
% 
% % Inizializzazioni
% cNonStatTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output degli intervalli
%                                       %(time frames) NON stazionari
% cMeanTacoTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output della media del taco
%                                        %nei vari time frames
% cStdTacoTimeFrames  = cell([1,nfasi]); %Inizializzo il cellarray di output della std del taco
%                                        %nei vari time frames
% cCImeanTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output degli intervalli di
%                                      %confidenza sulla media, calcolati per ogni fase della
%                                      %prova (un intervallo per ogni fase)
% cCIstdTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output degli intervalli di
%                                     %confidenza sulla dev.st., calcolati per ogni fase della
%                                     %prova (un intervallo per ogni fase)
% 
% % Ciclo su tutte le fasi della prova (length(timeFrames))
% for c = 1:length(timeFrames)
%     
%     % Sulla singola fase, verifico la stazionarieta' delle epoche di
%     % segnale considerate:
%     [rejFrames,vMean,vStd,ciMean,ciStd]=...
%         stationarityCheck(taco(:,2),taco(:,1),timeFrames{c},stationarityChkOptions);
%     
%     % Salvo gli estremi (tempo di inzio e fine, in ms) delle epoche NON
%     % stazionarie nell'apposito cell array di output e le elimino da
%     % "timeFrames" perche' da non considerare da qui in avanti:
%     cNonStatTimeFrames{c} = timeFrames{c}(rejFrames,:);
%     timeFrames{c} = timeFrames{c}(~rejFrames,:);
%     
%     % Salvo media e deviazione standard calcolate in ogni epoca in appositi
%     % cell array restituiti in output:
%     cMeanTacoTimeFrames{c} = vMean;
%     cStdTacoTimeFrames{c} = vStd;
%     
%     % Salvo gli intervalli di confidenza di media e deviazione standard,
%     % calcolati per la singola fase, in appositi cell array restituiti in
%     % output (utile per plot e verifiche):
%     cCImeanTimeFrames{c} = ciMean;
%     cCIstdTimeFrames{c} = ciStd;
%     
%     
% end


%--------------------------------------------------------------------------
%3. Stima dell'ordine ottimo per l'intero tacogramma passato in ingresso
%   (se parametro in input order NON e' scalare)

cOrdOptTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output dell'ordine
                                     %ottimo ottenuto per ciascun time-frame (ossia
                                     %per ciascuna finestra definita tramite
                                     %epochLen e overlap su tutte le fasi della
                                     %prova)
cPosOrdOptTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output del
                                        %minimo AIC scelto (il primo, il secondo,
                                        %il terzo, ecc.) combinando Aikake con 
                                        %test di bianchezza dei residui
%NB: Se "order" in input e' uno scalare, questi cell array resteranno vuoti
%    perche' non sara' stata effettuata alcuna selezione dell'ordine
%    ottimo.

if ~isscalar(order)
                                            
    % Ciclo su tutte le finestre individuate in timeFrames (cioe'
    % considerando sia epochLen che overlap), ottenendo gli ordini ottimi
    % per ciascuna epoca, all'interno di ciascuna fase della prova
    for c = 1:length(timeFrames)
        nepochs = size(timeFrames{c},1);
        vOrd = zeros(1,nepochs);
        vPosOrdOpt = zeros(1,nepochs);
        for ep = 1:nepochs
            %Identifico i campioni del taco di interesse per la singola epoca:
            p = taco(:,1)>=timeFrames{c}(ep,1) & taco(:,1)<=timeFrames{c}(ep,2);
            % Ciclo su tutti gli ordini possibili
            AICvalues = zeros(1,length(order));
            vAtest_h = zeros(1,length(order));
            vAtest_p = zeros(1,length(order));
            countOrd = 1;
            for ordAtt = order
                % Calcolo i coefficienti del modello AR di ordine p sulla
                % singola epoca di taco:
                armodel = ar(taco(p,2),ordAtt, 'yw') ;
                
                % Calcolo la cifra di Aikake per l'ordine considerato
                AICvalues(countOrd)=aic(armodel); %segmento � la 'finestra' di campioni che scelgo
                
                % Calcolo i residui (errore di predizione) con il modello
                % ottenuto:
                err = resid(taco(p,2),armodel);
                err = err.OutputData;
                
                % Verifico l'ipotesi di bianchezza dei residui con il test
                % di Anderson e ne memorizzo il risultato. 
                alphaAnders = residWhiteTestOpt.alpha;
                metAnders   = residWhiteTestOpt.method;
                [vAtest_h(countOrd),vAtest_p(countOrd)] = andersonWhiteTest(err,alphaAnders,metAnders);
                
                % Aggiorno il contatore del numero di ordini processati
                countOrd = countOrd+1;
            end
            
            % Ordino tutti le cifre di Aikake ottenute dal piu' piccolo al
            % piu' grande, conservandone le posizioni originali (per
            % ricordarmi a quale ordine si riferiscono). Riordino
            % coerentemente anche i risultati del test di Anderson e il
            % vettore degli ordini utilizzati:
            [~,P] = sort(AICvalues);
            vAtest_h = vAtest_h(P);
            vAtest_p = vAtest_p(P);
            orderSorted = order(P);
            
            % Trovo ordine ottimo per la singola epoca
            posOrdOpt = find(vAtest_h==false,1);
            
            % Se non viene trovato nessun ordine che soddisfa test di
            % Anderson, scelgo ordine associato a pValue maggiore (il piu'
            % vicino a soddisfare l'ipotesi di bianchezza dei residui) e
            % invio warning relativo:
            if isempty(posOrdOpt)
                warning('powerUnivariateAR:whitenessTestNotSatis',...
                        'Whiteness test not satisfied: returned optimal order with higher p-value');
                if strcmpi('adtest',residWhiteTestOpt.method)
                    
                    if max(vAtest_p) <= 4.9587 * 10^-5
                        % Condizione di sicurezza: se tutti i p-value sono
                        % pari o minori del valore minimo tabulato (4.9587*10^-6
                        % ma uso 10^-5 per sicurezza), allora non
                        % guardo piu' ai p-value, ma scelgo l'ordine
                        % massimo possibile:
                        [~,posOrdOpt] = max(orderSorted);
                    else
                        % Altrimenti agisco come spiegato sopra:
                        [~,posOrdOpt] = max(vAtest_p);
                    end
                    
                elseif strcmpi('savaresi',residWhiteTestOpt.method)
                    [~,posOrdOpt] = min(vAtest_p);
                end
            end
            
            % Memorizzo ordine ottimo per la singola epoca analizzata:
            vOrd(ep) = orderSorted(posOrdOpt);
            
            % Quale AIC minimo e' stato scelto per la singola epoca (il primo,
            % il secondo, il terzo, ecc.)?
            vPosOrdOpt(ep) = posOrdOpt;
            
        end
        
        % Salvo gli ordini ottimi trovati per tutte le epoche nella singola
        % fase della prova:
        cOrdOptTimeFrames{c} = vOrd;
        
        % Salvo la posizione del AIC minimo (il primo, il secondo, il
        % terzo, ecc.) scelto nella singola fase della prova considerando
        % anche l'esito del test di bianchezza:
        cPosOrdOptTimeFrames{c} = vPosOrdOpt;
        
    end
    
    % Calcolo l'ordine ottimo per il taco dato in ingresso:
    v = [];
    for c = 1:length(timeFrames)
        v = [v,cOrdOptTimeFrames{c}]; %#ok<AGROW>
    end
    % NB: con tante epoche a disposizione, le due opzioni sottostanti
    % dovrebbero produrre risultati confrontabili!
    %------------------------------------------
%     % OPZIONE A: Trovo l'ordine PIU' RICORRENTE tra quelli ottenuti per
%     %            tutte le epoche di segnale analizzate
%     catList = unique(v);
%     vCat = categorical(v,catList);
%     NCat = histcounts(vCat); %Numero di elementi per ogni "categoria" (cioe',
%                              %per ogni ordine effettivamente riconosciuto come
%                              %minimo per almeno un'epoca)
%     [~,p] = max(NCat);
%     optimalOrder = catList(p); %Ordine "ottimo" per il taco dato in input
    %------------------------------------------
    % OPZIONE B: Calcolo l'ordine MEDIANO tra quelli ottenuti per tutte le
    %            epoche di segnale analizzate
    optimalOrder = round(median(v));
    %------------------------------------------
    
    % Copio l'ordine ottimo scelto per il taco in input nella variabile
    % "order", per effettuare stima finale RSA con ordine ottimo ottenuto:
    order = optimalOrder;
end


%--------------------------------------------------------------------------
% 5. Calcolo la potenza del segnale nella banda HF considerata che coincide
%    in epoche di durata, in secondi, pari a "epochLen" e le medio tra
%    loro, restituendo un unico valore di RSA (con eventuale overlap tra le
%    epoche, se richiesto) per ogni fase della prova.
%    Effettuo anche test di bianchezza di Anderson per verificare che
%    ordine scelto sia effettivamente adatto alla maggior parte delle
%    epoche di segnale analizzate (utile anche per verificare l'efficacia
%    del metodo di esclusione degli spezzoni di segnale non stazionari) e
%    restituisco cell-array contenente H e p-value risultante dal test, per
%    ogni epoca considerata (una cella per ogni fase della prova, un
%    elemento del vettore per ogni epoca)

% Calcolo il valore medio di RSA per la fase considerata, utilizzando l'ordine "order"
% e considerando solo le epoche di segnale stazionario (cioe' quelle non
% scartate dal cell-array "timeFrames" al punto 2)
nfreqbands = size(freqbands,1);

BandPower = zeros(nfreqbands,nfasi); %Inizializzo la matrice di output
cBandPowerTimeFrames = cell([1,nfasi]); %Inizializzo il cellarray di output contenente le
                                        %potenze calcolate in ciascuna banda richiesta 
                                        %(ogni cella contiene una matrice che ha, per riga, 
                                        %una banda di frequenze considerata e, per colonna, 
                                        %un'epoca).
cPowerPolesTimeFrames = cell([1,nfasi]);%Inizializzo il cellarray di output contenente le
                                        %potenze associate a ciascun polo.
cFreqPolesTimeFrames = cell([1,nfasi]);%Inizializzo il cellarray di output contenente le
                                       %frequenze associate a ciascun polo.
cPSDtotTimeFrames = cell([1,nfasi]);  %Inizializzo il cellarray di output contenente la
                                      %PSD totale stimata sommando i contributi di tutti
                                      %i poli.
cAtestTimeFrames_h = cell([1,nfasi]); %Inizializzo il cellarray di output contenente la H del
                                      %test di bianchezza di Anderson:
                                      %H=0 --> i residui sono rumore bianco
                                      %H=1 --> i residui sono rumore colorato
cAtestTimeFrames_p = cell([1,nfasi]); %Inizializzo il cellarray di output contenente il
                                      %p-value del metodo di testing utilizzato 
                                      %(vedi funzione "andersonWhiteTest")

for c = 1:length(timeFrames)
    
    nepochs = size(timeFrames{c},1);
    vAtest_h = zeros(1,nepochs);
    vAtest_p = zeros(1,nepochs);
    mPot = zeros(nfreqbands,nepochs);
    mPowerPoles = NaN(order,nepochs); %NB: Inizializzo ad NaN perche', con stima 'onesided', il numero di poli restituiti
    mFreqPoles  = NaN(order,nepochs); %    dalla funzione "spectralDecompResiduals" puo' cambiare da epoca ad epoca (anche
                                      %    a parita' di ordine di modello AR!). Questo perche' il numero di poli
                                      %    complessi coniugati necessari per modellizzare il processo puo' cambiare tra i
                                      %    vari spezzoni di tacogramma.
    mPSDtot  = zeros(257,nepochs); %257 = numero di elementi dell'asse delle frequenze
                                   %      restituito dalla funzione "spectralDecompResiduals" 
                                   %      con stima 'onesided'.
    
    for ep = 1:nepochs
        %Identifico i campioni del taco di interesse per la singola epoca:
        p = taco(:,1)>=timeFrames{c}(ep,1) & taco(:,1)<=timeFrames{c}(ep,2);
        %Stimo i coefficienti del modello AR di ordine "order" fissato
        %(variabile scalare) e la VARIANZA DEL WHITE NOISE IN INGRESSO:
        [A,vep] = aryule(taco(p,2),order);
        %Calcolo i residui (errore di predizione) con il modello ottenuto:
        armodel = idpoly(A);
        err = resid(taco(p,2),armodel);
        err = err.OutputData;
        %Effettuo il test di bianchezza dei residui:
        alphaAnders = residWhiteTestOpt.alpha;
        metAnders   = residWhiteTestOpt.method;
        [vAtest_h(ep),vAtest_p(ep)] = andersonWhiteTest(err,alphaAnders,metAnders);
        
        %Calcolo la potenza riferita alla singola epoca (scegli uno):
        %----------------------------------------------
%         % METODO DI INTEGRAZIONE DELLO SPETTRO
%         nfft = 512; %Numero di punti su cui calcolare FFT (default = 512)
%         [psd,f] = pyulear(taco(p,2),order,nfft,fs,'onesided');
%         f_HF = f>=HFband(1) & f<=HFband(2);
%         vPot(ep) = bandpower(psd(f_HF),f(f_HF),'psd');
        %----------------------------------------------
        % METODO DEI RESIDUI (DECOMPOSIZIONE SPETTRALE)
        % Riferimento: Baselli 1997, script "sara.m" (che lo implementa
        %              sul modello AR monovariato), "mio_residui2.m" (by
        %              Giulia, utilizzato nella GUI; funzione adatta sia al
        %              modello mono che a quello multivariato).
        %
        [BandPowerEp,powerPoles,freqPoles,~,~,PSDtot,PSDfreqs] =...
           spectralDecompResiduals(1,A,vep,fs,freqbands,'onesided');
        mPot(:,ep) = BandPowerEp;
        mPowerPoles(1:length(powerPoles),ep) = powerPoles;
        mFreqPoles(1:length(freqPoles),ep) = freqPoles;
        mPSDtot(:,ep) = PSDtot;
        %----------------------------------------------
        
        
%         % VERIFICHE e PROVE (da commentare una volta testato il tutto):
%         [BandPowerEp,powerPoles,~,~,~,PSDtot,PSDfreqs] =...
%         spectralDecompResiduals(1,A,vep,fs,freqbands,'onesided');
%         [psd,f] = pyulear(taco(p,2),order,512,fs,'onesided');
%         %**************************
%         % Prova sul calcolo della potenza totale ottenuta con i due metodi
%         % (cioe' pyulear e metodo dei residui)
%         % RAZIONALE: La potenza totale (sia calcolata come area della PSD,
%         %            sia come somma della potenza associata ai singoli poli
%         %            (con il metodo dei residui), deve essere uguale alla
%         %            varianza del segnale in ingresso.
%         t = table(var(taco(p,2)),bandpower(psd,f,'psd'),sum(powerPoles),...
%             bandpower(PSDtot,PSDfreqs,'psd'),...
%             'VariableNames',{'Varianza','pyulear_potTot','somma_residui_mio','PSDtot_mio'});
%         disp(t);
%         %**************************
%         % Confronto potenza in banda ottenuta con i due metodi
%         % (cioe' integrazione dello spettro ["psd" ottenuta con pyulear e
%         % "PSDtot" dei residui sono del tutto equivalenti] vs potenza
%         % dei residui)
%         f_HF = PSDfreqs>=HFband(1) & PSDfreqs<=HFband(2);
%         t2 = table(bandpower(PSDtot(f_HF),PSDfreqs(f_HF),'psd'),...
%             BandPowerEp,...
%             'VariableNames',{'PotenzaInBanda_areaPSD','PotenzaInBanda_RESIDUI'});
%         disp(t2);
%         %**************************
%         % Plot diagnostico delle PSD ottenute
%         if ~exist('fPSD','var')
%             fPSD = figure;
%         end
%         figure(fPSD);
%         subplot(1,2,1);
%         plot(f,psd,PSDfreqs,PSDtot);
%         legend('pyulear','Residui mio');
%         title('Potenza NON normalizzata (CFR ESATTO)');
%         subplot(1,2,2);
%         plot(f,psd/bandpower(psd,f,'psd'),PSDfreqs,PSDtot/bandpower(PSDtot,PSDfreqs,'psd'));
%         legend('pyulear','Residui mio');
%         title('Potenza NORMALIZZATA (CFR MORFOLOGIA)');
%         %**************************
%         % Plot di confronto delle PSD ottenute con metodo parametrico
%         % (residui o pyulear, che in questo plot danno risultati IDENTICI)
%         % e NON parametrico (es. periodogramma di Welch)
%         if ~exist('fParam_vs_NonParam','var')
%             fParam_vs_NonParam = figure;
%         end
%         [psdW,f] = pwelch(taco(p,2),[],[],PSDfreqs,fs); % Metodo di confronto non parametrico
%         Y = fft(taco(p,2),512);
%         P2 = abs(Y); %double-sided spectrum (vedi "help fft"; non normalizzo per il numero di campioni perche' neanche gli altri metodi sono normalizzati)
%         P1 = P2(1:512/2+1); P1(2:end-1) = 2*P1(2:end-1); %single-sided spectrum (vedi "help fft"...)
%         P1plot = P1(1:512/2); % (vedi "help fft"...)
%         fFFT = 0:(fs/512):(fs/2-fs/512); % (vedi "help fft"...)
%         figure(fParam_vs_NonParam); plot(PSDfreqs,PSDtot,f,psdW,fFFT,P1plot);
%         legend('Residui o pyulear','Welch');
%         disp(['Fase = ',num2str(c),'   Epoca = ',num2str(ep)]);
%         %**************************
        
    end
    
    % Sistemo le matrici "mPowerPoles" e "mFreqPoles" in modo da
    % compattarle al massimo: rimuovo le righe delle due matrici costituite
    % da soli NaN
    s = sum(isnan(mPowerPoles),2);
    mPowerPoles(s==nepochs,:) = [];
    mFreqPoles(s==nepochs,:)  = [];
    
    % Salvo le variabili per singolo time frame negli appositi cell array:
    BandPower(:,c) = mean(mPot,2);
    cBandPowerTimeFrames{c} = mPot;
    cPowerPolesTimeFrames{c} = mPowerPoles;
    cFreqPolesTimeFrames{c} = mFreqPoles;
    cPSDtotTimeFrames{c} = mPSDtot;
    cAtestTimeFrames_h{c} = vAtest_h;
    cAtestTimeFrames_p{c} = vAtest_p;
end

%--------------------------------------------------------------------------
% 6. Preparo la struct "stats" da restituire in output:
% stats.cNonStatTimeFrames = cNonStatTimeFrames;
% stats.cStatTimeFrames = timeFrames; % E' giusto cosi'
% stats.cMeanTacoTimeFrames = cMeanTacoTimeFrames;
% stats.cStdTacoTimeFrames = cStdTacoTimeFrames;
% stats.cCImeanTimeFrames = cCImeanTimeFrames;
% stats.cCIstdTimeFrames = cCIstdTimeFrames;
stats.cOrdOptTimeFrames = cOrdOptTimeFrames;
stats.cPosOrdOptTimeFrames = cPosOrdOptTimeFrames;
stats.cBandPowerTimeFrames = cBandPowerTimeFrames;
stats.cPowerPolesTimeFrames = cPowerPolesTimeFrames;
stats.cFreqPolesTimeFrames = cFreqPolesTimeFrames;
stats.cPSDtotTimeFrames = cPSDtotTimeFrames;
stats.PSDfreqs = PSDfreqs;
stats.cAtestTimeFrames_h = cAtestTimeFrames_h;
stats.cAtestTimeFrames_p = cAtestTimeFrames_p;
stats.finalOrder = order;

end