A1 =  load('meancomplementariedad.txt'); 
A2 =  load('mincomplementariedad.txt'); 
A3 =  load('maxcomplementariedad.txt'); 
hold on
plot(A1(:,1),A1(:,2),'linewidth',3,'r--');
hold on
plot(A2(:,1),A2(:,2),'r--');
hold on
plot(A3(:,1),A3(:,2),'r--');
%xlabel(('Complementarity Threshold'),"fontsize",12)
%ylabel('Proportion Complementarity',"fontsize",12)
%print -color -F:12 complementarity.eps



