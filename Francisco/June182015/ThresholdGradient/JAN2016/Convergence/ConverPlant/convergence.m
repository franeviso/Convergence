OUT = zeros(100,10);
for t = 1:10;
PhenPlantR10;GenPlantR10;
G = 1 - G;
%add 0 to diag
G(logical(eye(size(G)))) = 0;
P(logical(eye(size(P)))) = 0;
%http://stackoverflow.com/questions/7781542/randperm-subset-of-random-m-by-n-matrix
m = 69;%size plant matrix
n = 69;

%# pick random rows
indX = randperm( size(G,1) );
indX = indX(1:m);

%# pick random columns
indY = indX;
indY = indY(1:n);

%# filter data
G = G(indX,indY);
P = P(indX,indY);

%remove sister species from GenAnimR
q = linspace(0,1,100);con = zeros(100,1);
for h = 1:100;
vector = 1:length(G);conver = 0;count = 0;count1=0;converA=0;converB=0;count2=0;converC=0;
for i = 1:length(G);sister = 0;nonsister=0;
i;
    %list of max gen sim
    maxvalue = [];
    [maxValue] = max(G(i,:));
    sister = find(G(i,:) == maxValue);
    %remove sister 
    nonsister = setdiff(vector,sister);
    nonsister(nonsister==i) = [];
    for j = 1:length(nonsister);
        P(i,:);
        G(i,:);
        P(i,sister);
        max(P(i,sister));
        P(i,nonsister(1,j));
        if P(i,nonsister(1,j)) > mean(P(i,sister));
        
           
            if P(i,nonsister(1,j)) > q(1,h);
               count = count + 1;
               converA = count;%/(length(nonsister));
            else
               count1 = count1 + 1;
               converB = count1;
            end
        else
          count2 = count2 + 1;
           converC = count2;
        end
        %conver
        %pause
    end
    %pheno i-nonsister > 0.95
end
con(h,1) = converA/(converA+converB+converC);
%fid = fopen('convergenceAll.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),con);fclose(fid); 
end%h
OUT(:,t) = con(:,1);
hold on
plot(q,con,'r')
end%t
fid = fopen('convergenceAll.txt','a');fprintf(fid, [repmat('% 6f',1,size(OUT,1)), '\n'],OUT);fclose(fid);
