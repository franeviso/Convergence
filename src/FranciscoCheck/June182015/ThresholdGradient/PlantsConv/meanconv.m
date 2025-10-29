A1 =  load('convergenciaPR1.txt'); 
A2 =  load('convergenciaPR2.txt'); 
A3 =  load('convergenciaPR3.txt'); 
A4 =  load('convergenciaPR4.txt'); 
A5 =  load('convergenciaPR5.txt'); 
A6 =  load('convergenciaPR6.txt'); 
A7 =  load('convergenciaPR7.txt'); 
A8 =  load('convergenciaPR8.txt'); 
A9 =  load('convergenciaPR9.txt'); 
A10 =  load('convergenciaPR10.txt'); 
A11 =  load('convergenciaPR11.txt'); 
A12 =  load('convergenciaPR12.txt'); 
A13 =  load('convergenciaPR13.txt'); 
A14 =  load('convergenciaPR14.txt'); 
A15 =  load('convergenciaPR15.txt'); 
Q = ((A1(:,2)+A2(:,2)+A3(:,2)+A4(:,2)+A5(:,2)+A6(:,2)+A7(:,2)+A8(:,2)+A9(:,2)+A10(:,2)+A11(:,2)+A12(:,2)+A13(:,2)+A14(:,2)+A15(:,2)))/15;
Q1 = (A1(:,1)+A2(:,1)+A3(:,1)+A4(:,1)+A5(:,1)+A6(:,1)+A7(:,1)+A8(:,1)+A9(:,1)+A10(:,1)+A11(:,1)+A12(:,1)+A13(:,1)+A14(:,1)+A15(:,1))/15;
h = plot(Q1,Q,'b');
hold on


Allrep;
for w = 1:19;
    MINI(w,1) = min(ALL(w,2:16));
    MAXI(w,1) = max(ALL(w,2:16));
end
plot(Q1,MINI,'b--');
hold on
plot(Q1,MAXI,'b--');
hold on


xlabel(('Convergence Threshold'),"fontsize",12)
ylabel('% Evolutionary Convergence Plants',"fontsize",12)
print -color -F:12 converplants.eps



