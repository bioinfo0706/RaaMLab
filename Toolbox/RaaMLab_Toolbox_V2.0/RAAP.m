function [AAP] = RAAP(Seq,reduce)
% Compute position-based features of reduced amino acids. 
% They describe the dispersion of the probability position distribution of 
% each reduced amino acid along the protein sequences using the coefficient
% of variation.
% Parameters of function
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,

% transform the protein sequence into reduced sequence% ÊäÈë²ÎÊý
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
% tranform the reduced sequence into the numerical sequences according to the
% interval distances of reduced amino acids
RSeq_p={1,G};
for si=1:G
    RSeq_p{si}=findstr(RSeq,AA(si));
end
% if loop
for si=1:G
    if length(RSeq_p{si})>1
        RSeq_p{si}=[RSeq_p{si},N+RSeq_p{si}(1)];
    end
end
% compute the position distribution of the reduced amino acids
mean=zeros(1,G);
add=zeros(1,G);
variance=zeros(1,G);
norm=zeros(1,G);
for si=1:G
    f=[];
    fre0=zeros(1,N);
    c={1,G};
    for sj=1:length(RSeq_p{si})-1
        c{si}(sj)=RSeq_p{si}(sj+1)-RSeq_p{si}(sj);           
        fre0(c{si}(sj))=fre0(c{si}(sj))+1;            
        f=fre0(:,1);                      
    end
    p=zeros(length(f));
    g=zeros(length(f));
    for st=1:length(f)
        p(st)=f(st)/length(c{si});      
        g(st)=st;
        mean(si)=mean(si)+g(st).*p(st);   
        add(si)=add(si)+(g(st)^2)*p(st);    
        variance(si)=add(si)-(mean(si)^2);  
        norm(si)=mean(si)/sqrt(variance(si));  
    end
    norm(isnan(norm))=0;
    AAP=norm;
end