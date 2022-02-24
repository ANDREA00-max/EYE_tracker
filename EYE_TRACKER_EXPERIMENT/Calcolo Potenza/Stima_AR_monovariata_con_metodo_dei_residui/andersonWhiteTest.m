function [h,p] = andersonWhiteTest(x,alpha,testMethod)

%%%% DESCRIZIONE:
% Implementazione test di bianchezza di Anderson, sia secondo la logica
% semplificata spiegata da Vercellis a lezione, sia utilizzando test di
% normalita' di Anderson-Darling disponibile su Matlab.
% Riferimenti: Savaresi, dispensa 07b, pp.8-12
%              "help resid"
%
%%%% USO:
% [h,p] = andersonWhiteTest(x,alpha,testMethod)
%
%%%% INPUTS:
% Solo il PRIMO E' OBBLIGATORIO, gli altri sono tutti OPZIONALI!
%
% x = vettore (riga o colonna) del segnale la cui bianchezza deve essere
%     testata.
%
% alpha = coefficiente di significativita' con cui si desidera condurre il
%         test di bianchezza.
%
% testMethod = metodo con cui si desidera condurre il test di aderenza della
%              ACF alla distribuzione normale teorica che ci si
%              aspetterebbe dalla ACF di un rumore bianco (escludendo
%              ACF(tau=0)), specificato come "savaresi" o "adtest" (caps
%              insensitive).
%
%%%% OUTPUTS:
%
% h = accettazione (h=0) o rifiuto (h=1) dell'ipotesi nulla di normalita'
%     della distribuzione della ACF del segnale "x" dato in input.
%     Corrisponde ad accettazione (h=0) o rifiuto (h=1) dell'ipotesi nulla
%     di bianchezza del segnale "x".
%
% p = a seconda del metodo, percentuale di valori della ACF di "x" che
%     stanno al di fuori dell'intervallo [-beta,+beta] (metodo "savaresi"),
%     oppure p-value del test di Anderson-Darling (metodo "adtest").
%
%%%%% INFORMAZIONI:
% Autore: Pierluigi Reali, PhD candidate
%         Dipartimento di Elettronica Informazione e Bioingegneria (DEIB)
%         Politecnico di Milano
% Anno: 2018
% Email: pierluigi.reali@polimi.it
%
%**************************************************************************

% Controllo parametri in input e inizializzazioni:
if nargin<2 || isempty(alpha)
    alpha = 0.05;
elseif alpha<=0 || alpha >=1
    error('andersonWhiteTest:alphaChk','Parametro alpha deve essere compreso tra 0 e 1');
end
if nargin<3 || isempty(testMethod)
    testMethod = 'savaresi';
else
    testMethod = lower(testMethod);
end

%--------------------------------------------------------------------------
% Numero di campioni del segnale in ingresso:
N = length(x);

%--------------------------------------------------------------------------
% Definisco il massimo lag (cioe', il massimo tau) da calcolare per la
% funzione di autocorrelazione (ACF).
% M=25 e' il numero che, di default, viene considerato dalla funzione
% matlab "resid" (System Identification Toolbox) per l'analisi dei residui.
% Mi sembra un numero un po' basso, pero', per effettuare un test ...
% Savaresi afferma che e' sufficiente considerare M<<N-1, mentre il test di
% Anderson-Darling, per dare risultati piu' robusti, richiede che il numero
% di campioni analizzati sia maggiore di 120 (in questo modo e' possibile
% utilizzare la "limiting distribution" definita da Anderson e Darling
% nell'articolo originale [vedi "help adtest"]). Percio' consideriamo (mia
% scelta motivata dalle ragioni precedenti) il numero di lag maggiore tra
% il 10% dei campioni di "x" e il valore costante 121.
% Riferimenti: "help resid"
%              "help adtest"
%              Savaresi, dispensa 07b, p. 12
M = max([121,round(0.10*N)]);

% Se il segnale X in input e' troppo corto rispetto al numero di lag di
% default, calcolero' la ACF per meno lags:
if M > N-1
    M = N-1;
end

%--------------------------------------------------------------------------
% Calcolo la funzione di auto-correlazione del segnale X (normalizzata):
[rho,tau] = xcorr(x,M,'coeff');

%--------------------------------------------------------------------------
% Calcolo l'indice che, dalla teoria (vedi Savaresi), dovrebbe essere
% distribuito come una normale N(0,1) (cioe' media=0 e dev.std=1):
rhoTest = rho(tau>=1)*sqrt(N);

%--------------------------------------------------------------------------
% Confronto la distribuzione di "rhoTest" con quella di una normale N(0,1).
% Posso usare il metodo suggerito da Savaresi nelle sue dispense (meno
% conservativo, piu' riluttante a considerare X un rumore bianco) oppure il
% test di normalita' di Anderson-Darling, con distribuzione di probabilita'
% specifica. Il punto di contatto teorico tra i due metodi e' spiegato
% all'inizio del case 'savaresi'.

switch (testMethod)
    case 'savaresi'
        %------------------------------------------------------------------
        % Definisco il valore di "beta", ossia il valore estremale della
        % distribuzione normale teorica N(0,1) al di la' del quale non ci
        % aspetteremmo di trovare più di alpha(%) valori.
        % Notiamo che il test e' svolto soltanto sulle code della
        % distribuzione, perche' e' li' che, evidentemente, ci si aspetta
        % la maggior deviazione dalla normalita' nel campo
        % dell'identificazione dei modelli: il test di normalita' di
        % Anderson-Darling ["help adtest"], infatti, e' sensibile
        % specialmente alle deviazioni sulle code (e, forse, e' da li' che
        % deriva il nome "test di bianchezza di Anderson" per il test
        % spiegato dal prof. Savaresi)
        mu = 0;    % Media della distribuzione normale teorica cercata
        sigma = 1; % Dev. standard        ""         ""           ""
        beta = norminv(1-(alpha/2),mu,sigma);
        
        %------------------------------------------------------------------
        % Contiamo il numero di punti di "rhoTest" che stanno al di fuori
        % dell'intervallo [-beta, +beta]:
        P = sum(rhoTest > beta | rhoTest < -beta);
        p = P/M; % ...in percentuale sul numero di lags calcolati.
        % NB: per definizione, si ha M=length(rhoTest)
        
        %------------------------------------------------------------------
        % Effettuiamo il test ed estraiamo H (rifiuto vs accettazione
        % ipotesi nulla)
        if p>alpha
            % Se la percentuale di valori che sta fuori da [-beta,+beta] e'
            % maggiore di quella che mi aspetterei (alpha), rifiuto l'ipotesi
            % nulla:
            h=true;
        else
            % ...altrimenti accetto l'ipotesi nulla:
            h=false;
        end
        
    case 'adtest'
        %------------------------------------------------------------------
        % Definisco l'oggetto "Probability Distribution" (vedi "help
        % makedist") associato alla normale teorica cercata N(0,1), cioe' a
        % media=0 e deviazione standard = 1
        pd = makedist('normal','mu',0,'sigma',1);
        
        %------------------------------------------------------------------
        % Eseguo il test di normalita' di Anderson-Darling (ipotesi nulla:
        % la distribuzione empirica testata NON differisce
        % significativamente dalla distribuzione teorica confrontata)
        [h,p]=adtest(rhoTest,'Distribution',pd,'Alpha',alpha);
        
    otherwise
        error('andersonWhiteTest:testMethodChk','Parametro testMethod non valido');
end


end