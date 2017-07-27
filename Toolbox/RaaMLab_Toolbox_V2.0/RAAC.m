function [SeqAAC,SeqDAAC,SeqTAAC,SeqFuAAC,SeqFiAAC] = RAAC(Seq,reduce)
% Calculate the content distribution of the reduced amino acids
% It contains reduced amino acid composition 
%             direduced-peptide composition, 
%             trireduced-peptide, 
%             tetrareduced-peptide composition
%             pentareduced-peptide compositionthe databases
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,

% tranform the protein sequence into reduced sequence
N=length(Seq);
AA='ARNDCQEGHILKMFPSTWYV';
G=max(reduce); 
RSeq=('');
for si=1:N
    for sj=1:20
        if Seq(si)==AA(sj)
            RSeq(si)=AA(reduce(sj));
        end
    end
end
% Compute the reuced amino acid composition
SeqAAC=zeros(1,G);
for si=1:G
    Aac=findstr(AA(si),RSeq);
    SeqAAC(1,si)=length(Aac)/N;
end
% Compute the direuced amino acid composition
DAA={G^2};
RAA=AA(1:G); 
n=1;
for si=1:G
    for sj=1:G
        DAA{n}=[RAA(si),RAA(sj)];
        n=n+1;
    end
end

SeqDAAC=zeros(1,G^2);
for si=1:G^2
    Aac=findstr(DAA{si},RSeq);
    SeqDAAC(1,si)=length(Aac)/(N-1);
end
% Compute the triuced amino acid composition
TAA={G^3};
RAA=AA(1:G); 
n=1;
for si=1:G
    for sj=1:G
        for sk=1:G
            TAA{n}=[RAA(si),RAA(sj),RAA(sk)];
            n=n+1;
        end
    end
end

SeqTAAC=zeros(1,G^3);
for si=1:G^3
    Aac=findstr(TAA{si},RSeq);
    SeqTAAC(1,si)=length(Aac)/(N-2);
end
% Compute the tetrareuced amino acid composition
FuAA={G^4};
RAA=AA(1:G);
n=1;
for si=1:G
    for sj=1:G
        for sk=1:G
            for st=1:G
                FuAA{n}=[RAA(si),RAA(sj),RAA(sk),RAA(st)];
                n=n+1;
            end
        end
    end
end

SeqFuAAC=zeros(1,G^4);
for si=1:G^4
    Aac=findstr(FuAA{si},RSeq);
    SeqFuAAC(1,si)=length(Aac)/(N-3);
end
% Compute the pentareuced amino acid composition
FiAA={G^5};
RAA=AA(1:G);
n=1;
for si=1:G
    for sj=1:G
        for sk=1:G
            for st=1:G
                for sr=1:G
                    FiAA{n}=[RAA(si),RAA(sj),RAA(sk),RAA(st),RAA(sr)];
                    n=n+1;
                end
            end
        end
    end
end

SeqFiAAC=zeros(1,G^5);
for si=1:G^5
    Aac=findstr(FiAA{si},RSeq);
    SeqFiAAC(1,si)=length(Aac)/(N-4);
end