function [Order_d,Quasi_order] = RSeqOrder(Seq,reduce,d,w)
% Compute two kinds of order-based features, one is sequence-order-coupling number, 
% and the other is quasi-sequence-order. They are computed based on the 
% Schneider–Wrede physicochemical distance matrix and Grantham chemical distance matrix.
% Order_d are results of sequence-order-coupling number, and
% Quasi_order shows results of quasi-sequence-order based on
% Schneider–Wrede physicochemical distance matrix and Grantham chemical distance matrix.
% Parameters of function
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% d-> the lag of the correlation
% w-> the weighting factor
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,


%transform the protein sequence into reduced sequence
matrix1=importdata('Grantham-distance.txt');
matrix2=importdata('Schneider-Wrede-distance.txt');
N=length(Seq); %求出序列长度
AA='ARNDCQEGHILKMFPSTWYV'; %给出氨基酸字母表
G=max(reduce); %求约化后字符的种类个数
RSeq=zeros(1,N);
RSeq2=zeros(1,N);
for si=1:N
    for sj=1:20
        if Seq(si)==AA(sj)
            RSeq(si)=AA(reduce(sj)); %序列转换成约化序列Seq==RSeq
            RSeq2(si)=reduce(sj);
        end
    end
end
%compute the reduced distance matrix based on the reduced amino acids
a={1,G};
for si=1:G
    a{si}=findstr(reduce,si);
end
aamatrix1={G,G};
aamatrix2={G,G};
AAmatrix1=zeros(G,G);
AAmatrix2=zeros(G,G);
for si=1:G
    for sj=1:G
        n=1;
        for st=1:length(a{si})
            for sk=1:length(a{sj})
                aa=[a{si}(st),a{sj}(sk)];
                aamatrix1{si,sj}(n)=matrix1(aa(1),aa(2));
                aamatrix2{si,sj}(n)=matrix2(aa(1),aa(2));
                n=n+1;
            end
            AAmatrix1(si,sj)=mean(aamatrix1{si,sj});
            AAmatrix2(si,sj)=mean(aamatrix2{si,sj});
        end
    end
end
% compute the sequence-order-coupling number
Order_d1=zeros(1,d);
Order_d2=zeros(1,d);
for si=1:d
     Order_d1(si)=0;
     Order_d2(si)=0;
    for sj=1:N-si
        Order_d1(si)=Order_d1(si)+(AAmatrix1(RSeq2(sj),RSeq2(sj+si)))^2;
        Order_d2(si)=Order_d2(si)+(AAmatrix2(RSeq2(sj),RSeq2(sj+si)))^2;
    end
    Order_d1(si)=Order_d1(si)/N-si;
    Order_d2(si)=Order_d2(si)/N-si;
end
Order_d=[Order_d1,Order_d1];
% compute the reduced amino acid composition
SeqAAC=zeros(1,G);
for si=1:G
    Aac=findstr(AA(si),RSeq);
    SeqAAC(1,si)=length(Aac)/N;
end
% compute the quasi-sequence-order
%
Quasi_order1=zeros(1,d);
Quasi_order2=zeros(1,d);
for si=1:G
    Quasi_order1(si)=SeqAAC(si)/(1+w*sum(Order_d1));
    Quasi_order2(si)=SeqAAC(si)/(1+w*sum(Order_d2));
end
%
for si=G+1:d
    Quasi_order1(si)=w*Order_d1(si)/(1+w*sum(Order_d1));
    Quasi_order2(si)=w*Order_d2(si)/(1+w*sum(Order_d2));
end
Quasi_order=[Quasi_order1,Quasi_order2];