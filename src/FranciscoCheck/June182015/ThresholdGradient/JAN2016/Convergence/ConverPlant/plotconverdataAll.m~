q = linspace(0,1,100);
A1 =  load('convergenceAll.txt'); 
%A2 =  load('mincomplementariedad.txt'); 
%A3 =  load('maxcomplementariedad.txt'); 
hold on
plot(q,mean(A1),'r');
hold on
An = sort(A1,'descend');

plot(q,An(10,:),'r--');
hold on
plot(q,An(90,:),'r--');

xlabel(('Convergence Threshold'),"fontsize",12)
ylabel('% Evolutionary Convergence Animals',"fontsize",12)
%print -color -F:12 complementarity.eps
