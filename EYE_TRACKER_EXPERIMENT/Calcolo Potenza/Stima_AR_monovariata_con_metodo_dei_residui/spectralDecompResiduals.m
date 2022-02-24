function [powerBands,powerPoles,freqPoles,poles,PSDpoles,PSDtot,PSDfreqs] = ...
          spectralDecompResiduals(B,A,vep,fs,freqBands,estimateT)

%%%%% USO:
% [powerBands,powerPoles,freqPoles,poles,PSDpoles,PSDtot,PSDfreqs] = ...
%             spectralDecompResiduals(B,A,vep,fs,freqBands,estimateT)
%
%%%%% DESCRIZIONE e NOTA TEORICA:
% NB: I "residui" di cui si parla qui non hanno nulla a che vedere con
%     l'errore di predizione. L'uso che se ne fa nel paper di Baselli,
%     intuitivamente, è abbastanza chiaro (e viene spiegato nel seguito di
%     questa descrizione), ma sembra che si parta da una conoscenza
%     pregressa di tale metodo (Zetterberg, 1979): Baselli, infatti, non ha
%     inventato il metodo dei residui ma ne propone l'applicazione sui suoi
%     modelli fisiologici.
%
% Funzione che effettua la decomposizione spettrale con metodo dei residui,
% a partire dai coefficienti del modello AR/ARMA precedentemente stimato
% sulla base dei dati (il "metodo dei residui" è, quindi, un metodo
% PARAMETRICO di stima spettrale).
% Si tratta di un metodo di "decomposizione spettrale" in quanto consente
% di separare l'effetto di ciascun polo sul calcolo della densita'
% spettrale di potenza (PSD). Viene definito "metodo dei residui"
% riferendosi al fatto che tale metodo richiede il calcolo della DISTANZA
% TRA TUTTI I POLI considerati nel modello. Nello specifico, viene
% calcolato il RECIPROCO della distanza tra i poli.
% Che senso ha calcolare il reciproco del "residuo" (ossia della DISTANZA)
% tra tutte le coppie di poli?
% Intuitivamente, più una coppia di poli è vicina, maggiore e' la potenza
% che troveremo nell'intorno delle frequenze caratterizzanti le due coppie
% di poli --> per questo il metodo dei residui usa il RECIPROCO DELLA
% DISTANZA tra le coppie di poli per calcolare la potenza loro associata:
% si prende in considerazione un polo alla volta e si calcola la
% produttoria della distanza tra questo polo e tutti gli altri; più la
% produttoria è piccola (poli molto vicini), maggiore sarà il suo reciproco
% e, quindi, la potenza associata al singolo polo.
%
% Provare per credere? Confronta il modulo della risposta in frequenza dei
% seguenti modelli:
% A = [1.0000   -3.1316    4.3358   -2.9600    0.8940]
% fvtool(1,A);
% A = [1.0000   -2.2552    2.8905   -2.1257    0.8940];
% fvtool(1,A);
% Dal diagramma poli/zeri, si puo' notare che l'unica differenza tra i due
% modelli e' data dall'avvicinamento (in termini di radianti) tra una
% coppia di poli complessi coniugati e l'altra. Dal diagramma del modulo,
% si nota che piu' i due poli sono "vicini" e maggiore e' la potenza
% complessiva nell'intorno di 0.2 di frequenza normalizzata (intervallo
% [0,1] tra [0 Hz,f_Nyquist]).
%
%%%%% INPUTS:
% Solo i primi 4 sono obbligatori, gli altri sono opzionali!
%
% B = vettore (riga o colonna) dei coefficienti al numeratore della FdT del modello
%     considerato (parte MA del modello) --> radici di B = zeri!
% A = vettore (riga o colonna) dei coefficienti al denominatore della FdT del modello
%     considerato (parte AR del modello) --> radici di A = poli!
% vep = varianza del rumore bianco in ingresso: varianza dello scarto tra il
%       valore stimato dal modello AR e il valore vero del campione
%       corrispondente.
%       Se il modello AR e' stato costruito con la funzione "aryule":
%               [A,vep] = ar(y,order,'yw');
%               
% fs = frequenza di campionamento del segnale utilizzato per stimare il
%      modello, sulla base della quale verranno calcolate le potenze per le
%      bande di interesse.
% freqBands = array delle bande di frequenza all'interno delle quali si vuole
%             calcolare la potenza associata ai poli (0, 1 , o piu') che
%             cadono all'interno di ciascuna banda.
%             Su ciascuna RIGA --> 1 banda di frequenze;
%             Su ciascuna delle 2 COLONNE --> estremi inferiore (1) e
%                                             superiore (2) di ciascuna
%                                             banda di interesse.
%             Di default, sono impostate le bande suggerite da Task Force
%             1996 per le bande LF ed HF della HRV (soggetti sani adulti).
% estimateT = Metodo di stima spettrale. Si puo' scegliere tra "onesided" e
%             "centered" (il significato di tali opzioni e' del tutto
%             equivalente all'uso che se ne fa nella funzione "pyulear"):
%             - onesided --> potenza e PSD vengono calcolate nel solo
%                            semiasse positivo delle frequenze.
%             - centered --> potenza e PSD vengono restituite con
%                            l'originale ripartizione tra semiasse
%                            negativo e positivo. ES: una singola coppia di
%                            poli complessi coniugati dara' origine ad una
%                            componente di potenza a frequenza positiva
%                            (associata al polo a parte immaginaria >0) e
%                            una a frequenza negativa (associata al polo a
%                            parte immaginaria <0) ).
%             Esso influenza:
%             - calcolo della POTENZA associata ai poli (powerPoles)
%               (ricavata dai residui) e, di conseguenza, il calcolo della
%               potenza in banda (powerBands).
%             - calcolo della PSD associata ai poli (PSDpoles, PSDtot).
%
%%%%% OUTPUTS:
% powerBands = potenza calcolata, con il metodo dei residui, nelle bande di
%              frequenza richieste.
% powerPoles = POTENZA associata ad ogni singolo polo del modello
%              (indipendentemente dalle "freqBands" date in input alla
%              funzione).
%              La somma di queste potenze, indipendentemente dal metodo di
%              stima scelto ("estimateT"), restituisce la varianza del
%              segnale originale.
% freqPoles = frequenza associata ad ogni polo (in Hertz), definita
%             rispetto alla frequenza di campionamento ("fs") data in
%             input.
% poles = radici di tutti i poli individuati nel modello considerato
%         (soluzione della eq.: denominatore_della_FdT = 0)
% PSDpoles = Contributo di ciascun polo alla DENSITA' SPETTRALE DI POTENZA
%            totale (decomposizione di "PSDtot" tra tutti i poli).
%            E' utile soprattutto per mostrare la ripartizione della
%            potenza tra i diversi poli:
%                  figure;
%                  plot(PSDfreqs,PSDtot,'k'); hold on;
%                  for k = 1:size(PSDpoles,2)
%                      area(PSDfreqs,PSDpoles(:,k),'FaceAlpha',0.3);
%                  end
%                  hold off;
%
% PSDtot = Densita' spettrale di potenza complessiva (del tutto equivalente
%          al risultato fornito da pyulear!), che coincide con la somma,
%          frequenza per frequenza, dei valori di "PSDtot".
% PSDfreqs = Asse delle frequenze sul quale sono stati valutati "PSDpoles"
%            e "PSDtot".
%
% RIFERIMENTI:
% Baselli 1997, script "sara.m" (che lo implementa sul modello AR
% monovariato), "mio_residui2.m" (by Giulia, utilizzato nella GUI; funzione
% adatta sia al modello mono che a quello multivariato).
%
%%%%% INFORMAZIONI:
% Autore: Pierluigi Reali, PhD candidate
%         Dipartimento di Elettronica Informazione e Bioingegneria (DEIB)
%         Politecnico di Milano
% Anno: 2018
% Email: pierluigi.reali@polimi.it
%
% Realizzato a partire da uno script di Giulia Tacchino (PhD.), integrato
% con possibilita' di scegliere tra stima "onesided" e "centered" e con
% possibilita' di estrarre piu' output (coerenti con il tipo di stima
% scelto).

%--------------------------------------------------------------------------
% 0. Inizializzazioni e verifiche sugli input

Ts = 1/fs; %Periodo di campionamento
freqBandsDef = [0.04,0.15; 0.15,0.4]; %Bande di frequenza di default
                                      %(LF e HF del segnale HRV, da Task
                                      %Force 1996)
estimateTDef = 'onesided';
estimateTValid = {'onesided','centered'};

if nargin < 5 || isempty(freqBands)
    freqBands = freqBandsDef;
end
if nargin < 6 || isempty(estimateT)
    estimateT = estimateTDef;
elseif sum( strcmpi(estimateT,estimateTValid) ) == 0
    error('spectralDecompResiduals:estimateTChkStr','Valore del parametro estimateT non valido');
else
    estimateT = lower(estimateT);
end

%---------
% CONTROLLO PROBABILMENTE DA ELIMINARE:
% % Per il momento, assumiamo che B (coefficienti del numeratore della FdT:
% % zeri) siano sempre pari a 1 --> modello AR monovariato: c'è solo la parte
% % autoregressiva.
% % Nel momento in cui vorremo aggiungere possibilità di calcolare la potenza
% % associata a modelli piu' complessi, allora dovremo rivedere il codice,
% % facendo sì che possa accettare in input anche B ~= 1.
% if B(1) ~= 1
%     error('spectralDecompResiduals:chkModelNum','Atteso primo coefficiente MA unitario (B(1)=1)');
% end
%---------

%--------------------------------------------------------------------------
% 1. Calcolo delle radici del polinomio al denominatore nel modello AR
%    e delle radici del polinomio al numeratore (se stiamo applicando la
%    decomposizione spettrale su un modello AR multivariato, assimilabile
%    ad un ARMA)
%          --> ottengo i POLI e gli ZERI del modello
poles = roots(A); %Poli del modello (vettore colonna)
if numel(B) ~=1
    % Se stiamo considerando modello ARMA (AR multivariato)...
    zeri  = roots(B); %Poli del modello (vettore colonna)
else
    % Se invece abbiamo a che fare con modello AR monovariato...
    zeri = [];
end

%--------------------------------------------------------------------------
% 2. Calcolo i "residui" (gamma) tra tutte le coppie di poli, considerando
%    anche eventuali zeri (nel caso di modello AR multivariato)
%    NB: Per comodita' e chiarezza, in questa sezione si e' trascurato
%        il contributo alla formula dei residui della varianza dell'errore
%        di predizione (vep): questo verra' aggiunto (cioe' moltiplicato)
%        nelle sezioni con il calcolo della potenza associata ai singoli
%        poli e della potenza spettrale totale e delle singole componenti.
%
% Riferimenti: formula in fondo a penultima pagina di Baselli 1997
%********************
% OPZIONE 1: METODO POCO EFFICIENTE, MA MOLTO CHIARO:
% Applichiamo direttamente la formula segnalata
gamma = zeros(length(poles),1); %Inizializzo il vettore dei residui:
                                %un residuo per ogni polo
for k = 1:length(poles) %Ciclo su tutti i poli e, per ognuno di essi (poles(k)),
                        %calcolo il residuo gamma(k)
    %-----------------
    % NUMERATORE DELLA FORMULA DEI RESIDUI (gamma(k)):
    if ~isempty(zeri)
        prod1 = 1; %Prima produttoria della formula
        prod2 = 1; %Seconda produttoria della formula
        for h = 1:length(zeri)
            prod1 = prod1*(poles(k) - zeri(h));
            prod2 = prod2*(poles(k)^-1 - zeri(h));
%             gamma_num = B(1) * prod1 * prod2; %Ciò che bisognerebbe
                                                %moltiplicare per le
                                                %produttorie e' il numeratore del
                                                %guadagno (che in questo
                                                %caso dovrebbe essere 1 ??)
            gamma_num = prod1 * prod2;
        end
    else
        gamma_num = complex(1);
    end
    %-----------------
    % DENOMINATORE DELLA FORMULA DEI RESIDUI (gamma(k)):
    prod1 =1; %Prima produttoria della formula
    prod2 =1; %Seconda produttoria della formula
    for h = 1:length(poles)
        if k~=h %La produttoria varrebbe 0 se svolgessimo anche: poles(k)-poles(h) ...
            prod1=prod1*(poles(k) - poles(h));
        end
        prod2=prod2*(poles(k)^-1 - poles(h));
    end
    gamma_den = poles(k) * prod1 * prod2; %Moltiplico il risultato delle due
                                          %produttorie tra loro e con poles(k), come indicato
                                          %nell'eq. in fondo alla penultima pagina di Baselli:
                                          % z * prod(z-ph) * prod(z^-1-ph);
    %-----------------
    % CALCOLO DEI RESIDUI:
    gamma(k) = gamma_num./gamma_den; %Calcolo i residui di ogni polo, mettendo insieme
                                     %numeratore e denominatore della formula 
end


%********************
% % OPZIONE 2: METODO PIU' ELEGANTE, ABBASTANZA CHIARO (autore: Pierluigi Reali)
% % Non offre alcun miglioramento, in termini di efficienza, rispetto al
% % metodo precedente (quindi lo tengo solo per riferimento futuro...)
% %
% % NB: L'operatore  .'  (punto+asterisco) esegue una trasposizione dei
% %     vettori complessi (scambia righe con colonne, come per i
% %     vettori reali); se utilizzassimo soltanto  '  (asterisco) con
% %     i vettori complessi, si otterrebbe, invece della trasposizione, il
% %     calcolo del complesso coniugato!
% Pk = repmat(poles.',order,1);
% Ph = repmat(poles,1,order);
% Pk_inv = Pk.^-1;
% PkPh = Pk-Ph;
% PkInvPh = Pk_inv-Ph;
% PkPh_I = PkPh+eye(order);
% vProd1 = prod(PkPh_I,1);
% vProd2 = prod(PkInvPh,1);
% res = poles.' .* vProd1 .* vProd2;
% res = res.^-1;
%********************
% % OPZIONE 3: METODO PIU' EFFICIENTE, MA POCO CHIARO (autore: Roberto Sassi)
% % Si tratta di una rielaborazione della formula di Baselli 1997, per
% % renderla piu' adatta alla risoluzione con il calcolatore). Ovviamente,
% % le tre opzioni forniscono risultati identici!
% % Riferimento: funzione "sara.m" di Roberto Sassi.
% gamma = prod((ones(order)+eye(order)-poles*((poles.^-1).')),1).* ...
%         prod((ones(order)-poles*(poles.')),1);
% gamma = gamma.^-1;
%********************


%--------------------------------------------------------------------------
% 3. Calcolo la potenza associata a ciascun polo, utilizzandone il residuo
%
% Serve anche per il calcolo della potenza nelle bande richieste. La
% somma della potenza associata a tutti i poli, sia con la stima
% "onesided" che "centered", restituisce la varianza (potenza totale)
% del segnale da cui e' stato stimato il modello AR --> esatto!!! E' logico
% e corrisponde a risultati forniti da "pyulear" e teoria (Baselli 1997).
% A seconda del tipo di stima spettrale richiesta, restituiamo frequenze
% negative ("centered") o soltanto positive ("onesided") e raddoppiamo la
% potenza associata ai poli complessi coniugati ("onesided"), oppure no
% ("centered").
% Qui si calcola anche dove si collocano i poli del modello, in termini di
% frequenza (Hz), sfruttando l'informazione fornita dalla frequenza di
% campionamento del segnale (la funzione "angle", di per sé, restituisce
% l'angolo di fase del polo in radianti [corrispondente all'angolo formato
% dal raggio passante per il polo, rispetto all'asse orizzontale, sulla
% circonferenza di raggio unitario]).
%
% Riferimenti: Penultima pagina di Baselli 1997, 2^ colonna

switch estimateT
    case 'onesided'
        % In una stima onesided non distinguiamo tra poli a parte
        % immaginaria positiva (angle(poles)>0, frequenza positiva) e a
        % parte immaginaria negativa (angle(poles)<0, frequenza negativa)
        freqPoles = abs( angle(poles)/(2*pi*Ts) );
        % Riordino poli e residui sulla base delle frequenze
        [freqPoles,I] = sort(freqPoles);
        poles = poles(I);
        gamma = gamma(I);
        % Cerco i poli complessi coniugati
        test = complex(real(poles),abs(imag(poles)));
        test = diff(test); %Gli 0 sono originati da poli con stessa parte
                           %reale e stesso modulo della parte immaginaria
        cc = find(test==complex(0));
        % Calcolo la potenza di tutti i poli
        % NB: Prendo la sola parte reale dei residui, come specificato in
        %     Baselli 1997, e non il modulo!
        powerPoles = real(gamma)*vep;
        % Raddoppio quella dei poli complessi coniugati ed eliminino le
        % componenti sul semiasse negativo delle frequenze (i poli puramente
        % reali, invece, verranno contati una volta sola). Raddoppio sia
        % perche' lo dice Baselli (1997), sia perche' in una stima "onesided",
        % siccome tutta la potenza del segnale deve concentrarsi sul
        % semiasse positivo delle frequenze, ha senso che la potenza dei
        % poli complessi coniugati venga riportata su quel semiasse
        % (necessario anche per garantire uguaglianza della potenza totale
        % [area sottesa allo spettro] nel caso di stima "onesided" e
        % "centered").
        powerPoles(cc) = powerPoles(cc)*2;
        powerPoles(cc+1) = [];
        freqPoles(cc+1) = [];

    case 'centered'
        % In una stima centered distinguiamo tra poli a parte immaginaria
        % positiva (angle(poles)>0, frequenza positiva) e a parte
        % immaginaria negativa (angle(poles)<0, frequenza negativa)
        freqPoles = angle(poles)/(2*pi*Ts);
        % Riordino poli e residui sulla base delle frequenze
        [freqPoles,I] = sort(freqPoles);
        poles = poles(I);
        gamma = gamma(I);
        % Calcolo la potenza di tutti i poli.
        % NB: Prendo la sola parte reale dei residui, come specificato in
        %     Baselli 1997, e non il modulo!
        powerPoles = real(gamma)*vep; 
end


%--------------------------------------------------------------------------
% 4. Calcolo la potenza nelle bande richieste
%    Se un polo "cade" nella banda attuale, allora la potenza del suo
%    residuo viene associata a tale banda.
powerBands = zeros(size(freqBands,1),1);

% Ciclo su tutte le bande di frequenza passate in input
for k = 1:size(freqBands,1)
    
    bandAtt = freqBands(k,:);
    
    % Determino quali poli cadono all'interno della banda considerata:
    % NB: NON considero il modulo della frequenza per coerenza con il tipo
    %     di stima richiesto in input ("onesided" o "centered")
    p = freqPoles>bandAtt(1) & freqPoles<=bandAtt(2);
    
    % Calcolo la potenza complessiva nella banda di interesse:
    % NB: "sum" di un vettore vuoto(nessun polo nel range di frequenze
    %     di interesse) restituisce 0 --> ok!!!
    powerBands(k) = sum( powerPoles(p) ); 
    
end


%--------------------------------------------------------------------------
% 5. Calcolo anche lo spettro di potenza associato a ciascun polo e lo
%    spettro totale che si otterrebbe considerando l'effetto complessivo di
%    tutti i poli (quest'ultimo e' lo stesso spettro che si otterrebbe con
%    funzione Matlab "pyulear", con opzioni "onesided" oppure "centered")

npoints = 512; %Numero di punti da considerare per lo spettro (= 1/risoluzione_in_frequenza)
               %Scegliere numero PARI per non avere incongruenze con pyulear 

% Inizializzazioni diverse in base a stima della PSD richiesta (scelte
% coerenti con quelle fatte dalla funzione "pyulear" per gli stessi valori
% di "nfft")
switch estimateT
    case 'onesided'
        %Spettro parziale relativo ad ogni singolo polo
        PSDpoles = zeros(npoints/2+1,length(poles));
        %Asse delle frequenze della PSD
        PSDfreqs = linspace(0,fs/2,npoints/2+1)';
    case 'centered'
        %Spettro parziale relativo ad ogni singolo polo
        PSDpoles = zeros(npoints,length(poles));
        %Asse delle frequenze della PSD
        PSDfreqs = linspace(-fs/2+1/npoints,fs/2,npoints)';
end

% Ciclo su tutti i poli
% NB: Le lettere "z", "S_segnato" e "S" usate nel ciclo sono coerenti con
%     quelle utilizzate in Baselli 1997 (penultima pagina)
for k = 1:length(poles)
    
    z = PSDfreqs*1i*2*pi*Ts;
    
    S_segnato = ( (gamma(k) * poles(k)) ./ (exp(-z)-poles(k)) ) +...
                gamma(k) + ...
                ( (gamma(k) * poles(k)) ./ (exp(z)-poles(k)) );
    
    S = S_segnato * Ts * vep; %I residui e i poli forniscono modulazione
                                 %dello spettro, mentre la varianza
                                 %dell'errore di predizione fornisce una base
                                 %("rescaling") dipendente dal segnale in ingresso
    PSDpoles(:,k) = S;
    
end

PSDpoles = real(PSDpoles); %Per il modulo dello spettro, devo tenere la sola
                           %parte reale (Baselli, 1997). Infatti, la parte
                           %immaginaria dello spettro dovrebbe essere
                           %legata alla fase...

% Nella stima "onesided", questa potenza raddoppia: l'area sottesa allo
% spettro (cioe' quella sottesa a "PSDtot"), infatti deve coincidere con la
% varianza del segnale originale (th. di Parseval), ma questa cosa deve
% valere qualunque sia la stima considerata e, quindi, sia tra 0 e pi che
% tra -pi e pi --> ho verificato: e' perfettamente coerente con
%                  comportamento di pyulear!
if strcmp(estimateT,'onesided')
    PSDpoles = PSDpoles*2;
end

PSDtot = sum(PSDpoles,2);  %Spettro complessivo (somma dei singoli effetti di tutti i poli)
                           %Risultato analogo a quello fornito da pyulear!
                           %Tutto corretto!    \(^_^)/
                           
end