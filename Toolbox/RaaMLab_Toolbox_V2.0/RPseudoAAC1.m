function [Psedo_AAC] = RPseudoAAC1(Index1,Index2,Index3,Seq,reduce,lambda,w)
% Calculate the type I pseudo-reduced amino acids of the protein sequences
% Index1 -> a kind of chosen physico-chemical properties of database1
% Index2 -> a kind of chosen physico-chemical properties of database1
% Index3 -> a kind of chosen physico-chemical properties of database1
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% lambda-> the lag of the RPseudoAAC;
% w-> the weighting factor
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,


% tansform the protein sequence into the reduced sequences
N=length(Seq); 
AA='ARNDCQEGHILKMFPSTWYV'; 
G=max(reduce); 
indextxt=importdata('AAindex.txt');
for si=1:length(indextxt.textdata)
     a=indextxt.textdata{si};
     b1=strcmp(Index1,a);
     b2=strcmp(Index2,a);
     b3=strcmp(Index3,a);
     if b1==1
         AAindex1=indextxt.data(si,:);
     end
     if b2==1
         AAindex2=indextxt.data(si,:);
     end
      if b3==1
         AAindex3=indextxt.data(si,:);
     end
end
ra={1,G};
rai1={1,G};
rai2={1,G};
rai3={1,G};
raii1=zeros(1,G);
raii2=zeros(1,G);
raii3=zeros(1,G);
for si=1:G
    ra{si}=findstr(reduce,si);
    for sj=1:length(ra{si})
        rai1{si}(sj)=AAindex1(ra{si}(sj));
        rai2{si}(sj)=AAindex2(ra{si}(sj));
        rai3{si}(sj)=AAindex3(ra{si}(sj));
    end
    raii1(si)=sum(rai1{si});
    raii2(si)=sum(rai2{si});
    raii3(si)=sum(rai3{si});
end
RAAindex1=zeros(1,20);
RAAindex2=zeros(1,20);
RAAindex3=zeros(1,20);
for si=1:G
    for sj=1:20
        if reduce(sj)==si
            RAAindex1(sj)=raii1(si);
            RAAindex2(sj)=raii2(si);
            RAAindex3(sj)=raii3(si);
        end
    end
end

RSeq=('');
for si=1:N
    for sj=1:20
        if Seq(si)==AA(sj)
            RSeq(si)=AA(reduce(sj)); 
        end
    end
end

RAA=AA(1:G); %the reduced amino acids
H1Seq=zeros(1,N);
H2Seq=zeros(1,N);
M1Seq=zeros(1,N);
for si=1:N
    for sj=1:G
        if Seq(si)==RAA(sj)
            H1Seq(si)=RAAindex1(sj);
            H2Seq(si)=RAAindex2(sj);
            M1Seq(si)=RAAindex3(sj);
        end
    end
end
% compute the reduced amino acid composition
SeqAAC=zeros(1,G);
for si=1:G
    Aac=findstr(AA(si),RSeq);
    SeqAAC(1,si)=length(Aac)/N;
end
% Compute Tk
Tk=zeros(1,lambda);
for si=1:lambda
    for sj=1:N-si
        Tk(si)=Tk(si)+((H1Seq(sj+si)-H1Seq(sj))^2+(H2Seq(sj+si)-H2Seq(sj))^2+(M1Seq(sj+si)-M1Seq(sj))^2)/3;
    end
    Tk(si)=Tk(si)/(N-si);
end

%compute the type I pseudo-reduced amino acids
Psedo_AAC=zeros(1,lambda);
for si=1:G
    Psedo_AAC(si)=SeqAAC(si)/(1+w*sum(Tk));
end
for si=G+1:G+lambda
    Psedo_AAC(si)=w*Tk(si-G)/(1+w*sum(Tk));
end