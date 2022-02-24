SNA_rat = [];
% clacolo SNA ratio a Low A

path=dir([pwd,'/Dati/S*']);

FRIED__HI = [];
FRIED__LOW = [];
FRI_VHF = [];
FRIED__ = [];

for j=1:length(path)
 
    %CARICO I DATI DEL FISSO 
    load([path(j).folder '/' path(j).name])

    for i = 1:3
    
   %calcolo delle bande di interesse
    
   %banda lf
  xLF = statsLM(i).PSDfreqs(5:21);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (5:21);

 ILF = trapz(xLF,yLF);
 
 %banda hf
 xHF = statsLM(i).PSDfreqs(21:53);
 yHF = statsLM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i,1)= ILF/IHF;
 SNA_HI(i,1) = IHF;
 SNA_LOW(i,1)= ILF;
 
 xLF = statsRM(i).PSDfreqs(5:21);
 yLF = statsRM(i).cPSDtotTimeFrames{1, 1} (5:21);

 ILF = trapz(xLF,yLF);
 
 
 xHF = statsRM(i).PSDfreqs(21:53);
 yHF = statsRM(i).cPSDtotTimeFrames{1, 1} (21:53);

 
 IHF = trapz(xHF,yHF);
 
 SNA_rat(i+3,1)= ILF/IHF; 
 SNA_HI(i+3,1) = IHF;
 SNA_LOW(i+3,1) = ILF;
end

% Arousal medio
 
for i = 4:6
    
   %calcolo delle bande di interesse
    
   %banda lf
  xLF = statsLM(i).PSDfreqs(5:21);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (5:21);

 ILF = trapz(xLF,yLF);
 
 %BANDA MA
 xHF = statsLM(i).PSDfreqs(21:53);
 yHF = statsLM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i-3,2)= ILF/IHF;
 SNA_HI(i-3,2) = IHF;
 SNA_LOW(i-3,2) = ILF;
 
 xLF = statsRM(i).PSDfreqs(5:21);
 yLF = statsRM(i).cPSDtotTimeFrames{1, 1} (5:21);

 ILF = trapz(xLF,yLF);
 
 
 xHF = statsRM(i).PSDfreqs(21:53);
 yHF = statsRM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i,2)= ILF/IHF;
 SNA_HI(i,2) = IHF;
 SNA_LOW(i,2) = ILF;
end
%BANDA HA

for i = 7:9
    
   %calcolo delle bande di interesse
    
   %banda lf
  xLF = statsLM(i).PSDfreqs(5:21);
  yLF = statsLM(i).cPSDtotTimeFrames{1, 1} (5:21);

 ILF = trapz(xLF,yLF);
 
   %banda hf
 xHF = statsLM(i).PSDfreqs(21:53);
 yHF = statsLM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i-6,3)= ILF/IHF;
 SNA_HI(i-6,3) = IHF;
 SNA_LOW(i-6,3) = ILF;
 
 xLF = statsRM(i).PSDfreqs(5:21);
 yLF = statsRM(i).cPSDtotTimeFrames{1, 1} (5:21);

 ILF = trapz(xLF,yLF);
 
 
 xHF = statsRM(i).PSDfreqs(21:53);
 yHF = statsRM(i).cPSDtotTimeFrames{1, 1} (21:53);

 IHF = trapz(xHF,yHF);
 
 SNA_rat(i-3,3)= ILF/IHF;
 SNA_HI(i-3,3) = IHF;
 SNA_LOW(i-3,3) = ILF;
end



 FRIED__ = [FRIED__;SNA_rat];
 FRIED__LOW = [FRIED__LOW,; SNA_LOW];
 FRIED__HI = [FRIED__HI; SNA_HI];
end


 
  
 
 
 
