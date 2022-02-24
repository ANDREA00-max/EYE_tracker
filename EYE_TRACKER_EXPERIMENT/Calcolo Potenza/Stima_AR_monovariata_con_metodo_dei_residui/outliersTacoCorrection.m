function [tacoCorr,vOut] = outliersTacoCorrection(taco,nSD)

%%%%% INPUTS
% taco = vettore(riga o colonna)/matrice del tacogramma che contiene:
%          - se VETTORE: la sola distanza temporale tra picchi R
%                        consecutivi (SENZA l'istante di rilevazione di
%                        ciascun picco R).
%          - se MATRICE (n x 2): 1^ colonna: istante di rilevazione di
%                                ciascun picco R;
%                                2^ colonna: distanza temporale tra picchi
%                                R consecuvitvi.
% nSD  = distanza dalla media in termini di deviazione standard che si
%        vuole utilizzare come criterio per l'identificazione degli
%        outliers. Kemper 2007 suggerisce di utilizzare nSD=5, ma si
%        riferisce a nati prematuri, mentre Tacchino 2011 sceglie nSD=4 per
%        adulti.
%
%%%%% OUTPUTS
% tacoCorr = vettore/matrice (coerentemente con l'input "taco") delle
%            distanze temporali tra picchi R consecutivi corrette con il
%            metodo di sostituzione degli outliers implementato.
% vOut     = righe della matrice/elementi del vettore del tacogramma di
%            partenza in cui sono stati individuati outliers.
%
%%%%% RIFERIMENTI:
% Kemper 2007, Tacchino 2011, Barnaby 2002
%
%%%%% INFORMAZIONI:
% Autore: Pierluigi Reali, PhD candidate
%         Dipartimento di Elettronica Informazione e Bioingegneria (DEIB)
%         Politecnico di Milano
% Anno: 2018
% Email: pierluigi.reali@polimi.it

%**************************************************************************
% Controlli e Inizializzazioni
flag_t = false; % inizializzazione flag di trasposizione dell'input
flag_v = false; % inizializzazione flag vettore in input

if isvector(taco)
    % Se "taco" e' un vettore, assumero' campionamento a frequenza
    % costante, interpolero' gli outliers in tal senso (metodo "fillmissing") e
    % restituiro' un vettore tacoCorr (serie numerica) in uscita.
    flag_v = true;
    if isrow(taco)
        taco = taco';
        flag_t = true;
    end
    
elseif ismatrix(taco) && size(taco,2) == 2 %ismatrix testa se "taco" e' matrice 2D
    % Se "taco" e' una matrice, il campionamento puo' anche essere a
    % frequenza variabile: l'istante temporale a cui si riferisce ciascun
    % campione del segnale HRV e' contenuto nella 1^ colonna di "taco",
    % mentre il valore di distanza RR assunto dal segnale in quell'istante
    % e' contenuto nella 2^ colonna della corrispondente riga. Per ottenere
    % un risultato piu' accurato, interpolero' gli outliers basandomi sulla
    % scala temporale di partenza (metodo "interp1") e e restituiro' una
    % matrice tacoCorr (serie temporale) in uscita.
    
    vtime = taco(:,1);
    taco  = taco(:,2);
else
    error('outliersTacoCorrection:invalidTaco','Dimensione del input taco non valida.');
end


N = length(taco);
wlen = 100;  % Dimensione della finestra mobile (sia Kemper 2007 che
             % Tacchino 2011 suggeriscono wlen=100)
percent = 20; % Per come e' implementata la funzione trimmean (vedi help trimmean),
              % esempio:
              % 20 --> calcolo della media con esclusione del 10% dei
              %        campioni superiori (sopra il 90° percentile) e del
              %        10% dei campioni inferiori (sotto il 10° percentile)
              %        del dataset.
vOut = [];   % Inizializzo il vettore degli outliers da restituire in uscita.

%**************************************************************************
% Identificazione e sostituzione outliers sui primi 100 campioni del
% tacogramma (inizializzazione algoritmo di identificazione outliers)

% 0. Prendo lo spezzone di tacogramma da considerare per l'inizializzazione
%    dell'algoritmo di identificazione outliers
cEndInit = min([N,wlen]); %Per gestire correttamente i tacogrammi che hanno
                          %numero di battiti <= 100
tacoInit = taco(1:cEndInit);

% 1. Calcolo la media troncata sull'intero tacogramma
m = trimmean(taco,percent);

% 2. Identifico gli outliers sui primi 100 campioni con il criterio
%    suggerito

% vOutInit = find( tacoInit> m+m*0.67 | tacoInit< m-m*0.67); %Tacchino 2011
vOutInit = find( tacoInit> m*2 | tacoInit< m*0.5); %Kemper 2007

% 3. Inizializzo il tacogramma corretto a quello originale
tacoCorr = tacoInit;

% 4. Se sono stati individuati outliers in questa prima fase
%    (inizializzazione algoritmo), correggo tali intervalli RR
%    sostituendovi il valore trovato con la media troncata
tacoCorr(vOutInit) = m;

%**************************************************************************
% Identificazione e sostituzione outliers sui successivi campioni del
% tacogramma (algoritmo vero e proprio per identificazione outliers)

% Eseguo questa parte dell'algoritmo solo se il numero di campioni
% complessivi del tacogramma e' superiore a 100
if N > wlen
    
    t  = wlen+1; %Primo campione del taco da controllare:
                 %comincio dal campione successivo a quello di fine inizializzazione
    
    % IDENTIFICO gli outliers sul tacogramma
    while t <= N
        
        % 1. Calcolo media e deviazione standard degli intervalli RR
        %    sull'ultima porzione di tacogramma gia' corretta
        ti = t-wlen; %Campione iniziale della finestra mobile
        tf = t-1;    %Campione finale della finestra mobile
        
        mWin  = nanmean(tacoCorr(ti:tf));
        sdWin = nanstd(tacoCorr(ti:tf));
        
        % 2. Definisco il range per individuare gli outliers
        th = [mWin-nSD*sdWin, mWin+nSD*sdWin];
        
        % 3. Verifico se il campione attuale (t) e' un outlier oppure no e,
        %    nel caso, lo contrassegno come NaN:
        if taco(t) < th(1) || taco(t) > th(2)
            vOut = [vOut;t]; %#ok<AGROW>
            tacoCorr(t) = NaN;
        else
            tacoCorr(t) = taco(t);
        end
        
        % 4. Passo al campione successivo
        t = t+1;
        
    end
    
    % 6. SOSTITUISCO gli intervalli RR che sono outliers tramite interpolazione
    %    spline dei campioni accettati vicini, come suggerito in Kemper 2007 e
    %    Barnaby 2002
    if flag_v == true
        % Se taco in ingresso e' un vettore:
        tacoCorr = fillmissing(tacoCorr,'spline');
    else
        % Se taco in ingresso e' una matrice, con una specifica scala dei tempi:
        tacoCorrNoNaN = tacoCorr(~isnan(tacoCorr));
        vtimeNoNaN    = vtime(~isnan(tacoCorr));
        tacoCorr = interp1(vtimeNoNaN,tacoCorrNoNaN,vtime,'spline');
        tacoCorr = [vtime,tacoCorr];
    end
    
    
end

%**************************************************************************
% Memorizzo nel vettore in uscita vOut tutti gli outliers trovati
vOut = [vOutInit; vOut];

% Sistemo l'output sulla base della dimensionalita' del tacogramma in input
if flag_t ==true
    tacoCorr=tacoCorr';
end


end
