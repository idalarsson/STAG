
function [X,design] = datatobarcodedata(expMat)
%INPUT
%-----
%expMat: Should be table with cells as rows and columns are
%[genes,barcode,day,treatment,state]. Important that columns are arranged in
%correct order.
% 
%Formatting of expression matrix: 
%Barcode - 
%Days - 
%Treatment - Value from 1 and up. Control should always be 1. 
%
%OUTPUT
%-----
%X: a matrix of states x barcodes 
%design: a matrix specifying the design of the experiment
%These are to be used as input for netfrombarcode.m

[c,d]=size(expMat);
[state,~,s]=unique(expMat(:,d)); 
[treatment,~,tt]=unique(expMat(:,d-1));
[tp,~,t]=unique(expMat(:,d-2));
[bctag,~,bc]=unique(expMat(:,d-3)); 
bctag=table2cell(bctag);
x=size(state,1);
y=max(tt)*max(bc)*max(t);
X=zeros(x,y);

for i=1:max(tt)
    for j=1:max(t)
        for k=1:max(bc)
            for l=1:max(s)
                m=((i-1)*(max(bc)*max(t)))+((j-1)*max(bc))+k;
                f=intersect(intersect(intersect(find(bc==k),find(s==l)),find(t==j)),find(tt==i)); %add one intersect for treatments
                X(l,m)=length(f);
            end
        end
    end
end

design=zeros(2,y); 

for i=1:max(tt) 
    %repmat(i,1,max(bc)*max(t))
    design(1, max(bc)*max(t)*(i-1)+1 : max(bc)*max(t)*(i-1) + max(bc)*max(t))=repmat(i,1,max(bc)*max(t));
    for j=1:max(t)
        design(2 , max(bc)*(j-1)+1 +(max(bc)*max(t)*(i-1)) :  max(bc)*j+(max(bc)*max(t)*(i-1)))=repmat(j,1,max(bc));
    end
end

end