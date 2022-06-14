data = rand(5,5);
m = 3;
n = 3;

%# pick random rows
indX = randperm( size(data,1) )
indX = indX(1:m)

%# pick random columns
indY = randperm( size(data,2) )
indY = indY(1:n)

%# filter data
data2 = data(indX,indY)
