function [Composition,Transition,Distribution] = RCTD(Seq,reduce)
% Calculate the composition,transition,distribution reduced amino acids
%            composition: reduced amino acid composition
%            transition: A transition from a reduced amino acid  to another
%            is the percent frequency with which the reduced amino acid  is followed by the reduced amino acid  or A is followed by I in the reduced sequence. 
%            distribution: There are five ¡°distribution¡± descriptors for
%            each attribute and they are the position percents in the whole sequence for the first residue, 25% residues, 50% residues, 
%           75% residues and 100% residues , respectively, for a specified
%           encoded class.
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,

% tranform the protein sequence into reduced sequence
N=length(Seq); 
AA='ARNDCQEGHILKMFPSTWYV'; 
G=max(reduce);
RSeq=('');
RAA=AA(1:G); 
for si=1:N
    for sj=1:20
        if Seq(si)==AA(sj)
            RSeq(si)=RAA(reduce(sj)); 
        end
    end
end
% Compute the composition
Composition=zeros(1,G);
for si=1:G
    Aa=findstr(RAA(si),RSeq);
    Composition(1,si)=length(Aa)/N; 
end
%Compute the translation
RAA=AA(1:G); 
n=1;
aa={''};
for si=1:G
    for sj=si+1:G
        aa{1,n}=[RAA(si),RAA(sj)];
        aa{2,n}=[aa{1,n}(2),aa{1,n}(1)];
        n=n+1;
    end
end
Transition=zeros(1,size(aa,2));
for si=1:size(aa,2)
    t1=length(findstr(aa{1,si},RSeq));
    t2=length(findstr(aa{2,si},RSeq));
    Transition(si)=(t1+t2)/(N-1);
end
%Compute the distribution
RAA=AA(1:G); 
Aac={1,G};
for si=1:G
    Aac{si}=findstr(RAA(si),RSeq);
end
D={G};
for si=1:G
    if isempty(Aac{si})
        D{si}(1:5)=0;
%    elseif Aac{si}(1)==1
%        D{si}(1)=0;
    else
        D{si}(1)=Aac{si}(1)/N;
        D{si}(2)=Aac{si}(round(length(Aac{si})*0.25))/N;
        D{si}(3)=Aac{si}(round(length(Aac{si})*0.50))/N;
        D{si}(4)=Aac{si}(round(length(Aac{si})*0.75))/N;
        D{si}(5)=Aac{si}(length(Aac{si}))/N;
    end
end
Distribution=cell2mat(D);
%%