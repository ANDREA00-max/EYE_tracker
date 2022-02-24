SNA_rat = [];
% clacolo SNA ratio a Low A

path=dir([pwd,'/Dati/S*']);


FRIED_SNA = [];

for j=1:length(path)
     
    %CARICO I DATI DEL FISSO 
    load([path(j).folder '/' path(j).name '/' path(j).name '_arousal_con_PSD.mat'])
   

    for i = 1:3
    
   %calcolo delle bande di interesse
    
   %banda lf
  xLF = statsLM(i).PSDfreqs(1:21);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (1:21);

 ILF = trapz(xLF,yLF);
 
 %banda hf
 xHF = statsLM(i).PSDfreqs(21:53);
 yHF = statsLM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i,1)= ILF/IHF;
 
 xLF = statsRM(i).PSDfreqs(1:21);
 yLF = statsRM(i).cPSDtotTimeFrames{1, 1} (1:21);

 ILF = trapz(xLF,yLF);
 
 
 xHF = statsRM(i).PSDfreqs(21:53);
 yHF = statsRM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i+3,1)= ILF/IHF;
end

% Arousal medio
 
for i = 4:6
    
   %calcolo delle bande di interesse
    
   %banda lf
  xLF = statsLM(i).PSDfreqs(1:21);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (1:21);

 ILF = trapz(xLF,yLF);
 
 %BANDA MA
 xHF = statsLM(i).PSDfreqs(21:53);
 yHF = statsLM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i-3,2)= ILF/IHF;
 
 xLF = statsRM(i).PSDfreqs(1:21);
 yLF = statsRM(i).cPSDtotTimeFrames{1, 1} (1:21);

 ILF = trapz(xLF,yLF);
 
 
 xHF = statsRM(i).PSDfreqs(21:53);
 yHF = statsRM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i,2)= ILF/IHF;
end
%BANDA HA

for i = 7:9
    
   %calcolo delle bande di interesse
    
   %banda lf
  xLF = statsLM(i).PSDfreqs(1:21);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (1:21);

 ILF = trapz(xLF,yLF);
 
   %banda hf
 xHF = statsLM(i).PSDfreqs(21:53);
 yHF = statsLM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i-6,3)= ILF/IHF;
 
 xLF = statsRM(i).PSDfreqs(1:21);
 yLF = statsRM(i).cPSDtotTimeFrames{1, 1} (1:21);

 ILF = trapz(xLF,yLF);
 
 
 xHF = statsRM(i).PSDfreqs(21:53);
 yHF = statsRM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i-3,3)= ILF/IHF;
end


 FRIED_SNA = [FRIED_SNA;SNA_rat];

 
  
 
 
 