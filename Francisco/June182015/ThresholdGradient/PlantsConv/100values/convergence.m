PhenPlantsR15;GenPlantsR15;
%add 0 to diag
G(logical(eye(size(G)))) = 0;
%remove sister species from GenAnimR
for h = 1:100;
%q = [0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98];
q = linspace(0,1,100);
vector = 1:length(G);converA = 0;converB=0;count = 0;count1=0;
for i = 1:length(G);
i;
    %list of max gen sim
    maxvalue = [];
    [maxValue] = max(G(i,:));
    sister = find(G(i,:) == maxValue);
    %remove sister 
    nonsister = setdiff(vector,sister);
    nonsister(nonsister==i) = [];
    %compare covergence with non-sister sp using a threshold
    %pheno simil i-nonsister > pheno i-sister
    for j = 1:length(nonsister);
        P(i,:);
        G(i,:);
        P(i,sister);
        max(P(i,sister));
        P(i,nonsister(1,j));
        if P(i,nonsister(1,j)) > max(P(i,sister));
           
            if P(i,nonsister(1,j)) > q(1,h);
               count = count + 1;
               converA = count;%/(length(nonsister));
            else
               count1 = count1 + 1;
               converB = count1;
            end
        end
        %conver
        %pause
    end
    %pheno i-nonsister > 0.95
end
con = converA/(converA+converB);
fid = fopen('convergenciaPR15q100.txt','a');fprintf(fid,'%6.8f %6.8f\n',q(1,h),con);fclose(fid); 
%fid = fopen('OutputAll1000.txt','a');fprintf(fid, [repmat('% 6f',1,size(abup1,2)), '\n'],abup1);
end

