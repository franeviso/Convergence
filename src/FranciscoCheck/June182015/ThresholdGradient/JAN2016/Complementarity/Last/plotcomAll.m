q = linspace(0,1,100);
A1 =  load('complementarityAll.txt'); 
%A2 =  load('mincomplementariedad.txt'); 
%A3 =  load('maxcomplementariedad.txt'); 
hold on
plot(q,mean(A1),'r','LineWidth',1);
hold on
An = sort(A1,'descend');

plot(q,An(1,:),'r--');
hold on
plot(q,An(9,:),'r--');

xlabel(('Complementarity Threshold'),"fontsize",12)
ylabel('% Complementarity',"fontsize",12)
%print -color -F:12 complementarity.eps
