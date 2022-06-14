A1 =  load('convergenciaAR1.txt'); 
A2 =  load('convergenciaAR2.txt'); 
%A3 =  load('convergenciaAR3.txt'); 
%A4 =  load('convergenciaAR4.txt'); 
%A5 =  load('convergenciaAR5.txt'); 
%A6 =  load('convergenciaAR6.txt'); 
%A7 =  load('convergenciaAR7.txt'); 
%A8 =  load('convergenciaAR8.txt'); 
%A9 =  load('convergenciaAR9.txt'); 
%A10 =  load('convergenciaAR10.txt'); 
%A11 =  load('convergenciaAR11.txt'); 
%A12 =  load('convergenciaAR12.txt'); 
%A13 =  load('convergenciaAR13.txt'); 
%A14 =  load('convergenciaAR14.txt'); 
%A15 =  load('convergenciaAR15.txt'); 
Q = (A1(:,2)+A2(:,2))/2;
Q1 = (A1(:,1)+A2(:,1))/2;
%Q = ((A1(:,2)+A2(:,2)+A3(:,2)+A4(:,2)+A5(:,2)+A6(:,2)+A7(:,2)+A8(:,2)+A9(:,2)+A10(:,2)+A11(:,2)+A12(:,2)+A13(:,2)+A14(:,2)+A15(:,2)))/15;
%Q1 = (A1(:,1)+A2(:,1)+A3(:,1)+A4(:,1)+A5(:,1)+A6(:,1)+A7(:,1)+A8(:,1)+A9(:,1)+A10(:,1)+A11(:,1)+A12(:,1)+A13(:,1)+A14(:,1)+A15(:,1))/15;
h = plot(Q1,Q,'r');
hold on


%Allrep;
%for w = 1:2;
%    MINI(w,1) = min(ALL(w,2:16));
%    MAXI(w,1) = max(ALL(w,2:16));
%end
plot(Q1,A1(:,2),'r--');
hold on
plot(Q1,A2(:,2),'r--');
hold on


xlabel(('Convergence Threshold'),"fontsize",16)
ylabel('% Evolutionary Convergence Plants',"fontsize",16)
print -color -F:12 converplants.eps



