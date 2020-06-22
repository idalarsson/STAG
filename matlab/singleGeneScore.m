
function [states] = singleGeneScore(expM,gene,method)
%gene='GDF15';
tidx = find(string(expM.Properties.VariableNames) == string(gene))
[a,~]=size(expM);
states=zeros(a,1);

if strcmp(method,'mean')
    c=mean(table2array(expM(:,tidx)));
    for i=1:a
        if (table2array(expM(i,tidx)) > c)
            states(i,1)=2;
        else
            states(i,1)=1;
        end
    end
elseif strcmp(method,'quantile')
    c=quantile(table2array(expM(:,tidx)),[0.25 0.5 0.75]);
    for i=1:a
        if (table2array(expM(i,tidx)) < c(1))
            states(i,1)=1;
        elseif (table2array(expM(i,tidx)) > c(3))
            states(i,1)=4;
        elseif (c(1)<table2array(expM(i,tidx))) && (table2array(expM(i,tidx)) < c(2))
            states(i,1)=2;
        else
            states(i,1)=3;
            
        end
    end
end

states=array2table(states);
end
