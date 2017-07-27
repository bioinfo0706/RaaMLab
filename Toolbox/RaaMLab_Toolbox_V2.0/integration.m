function [] = integration(filename,Index,para_dis,para_meth,maxgroup)
%%
% Calculate the content distribution of the reduced amino acids
% It contains reduced amino acid composition 
%             direduced-peptide composition, 
%             trireduced-peptide, 
%             tetrareduced-peptide composition
%             pentareduced-peptide compositionthe databases
% 'filename' is a name of a multiple fasta file
% Index-> the index of databases in RaaMlab toolbox
% para_dis-> computes the distance of amino acids based on the index of databases
% para_meth-> creates a hierarchical cluster tree, using the single
% maxgroup_> is the size of reduced amino acid set. Maxgroup can be from 1 to 20.
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,
%% load dataset
% reads set of sequence from multyple fasta format
% 'filename' is a name of a multiple fasta file
[maxLen,nSeq,Seq,SeqName,SeqLength]=readfasta(filename);
%% AAreduce Parameters
% Index-> the index of databases in RaaMlab toolbox
% para_dis-> computes the distance of amino acids based on the index of databases
% para_meth-> creates a hierarchical cluster tree, using the single
% maxgroup_> is the size of reduced amino acid set. Maxgroup can be from 1 to 20.
reduce=AAreduce(Index,para_dis,para_meth,maxgroup);
%% K-mer
% Parameters of function
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids
[aac,daac,taac,fuaac,fiaac]=RAAC(Seq,reduce)
Output('RAAC_aac',aac);
Output('RAAC_daac',daac);
Output('RAAC_taac',taac);
Output('RAAC_fuaac',fuaac);
Output('RAAC_fiaac',fiaac);
%% Composition, transition and distribution
% Parameters of function
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids
[co,tr,di]=RCTD(Seq,reduce);
Output('RCTD_co',co);
Output('RCTD_tr',tr);
Output('RCTD_di',di);
%% pseudo-reduced amino acid composition (PRseAAC)
% PRseAAC Parameters
% Calculate the type I pseudo-reduced amino acids of the protein sequences
% Index1 -> a kind of chosen physico-chemical properties of database1
% Index2 -> a kind of chosen physico-chemical properties of database1
% Index3 -> a kind of chosen physico-chemical properties of database1
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% lambda-> the lag of the RPseudoAAC;
% w-> the weighting factor
Index1='ARGP820103';
Index2='BEGF750102';
Index3='BHAR880101';
lambda=10;
w=0.03;
pseudoaac1=RPseudoAAC1(Index1,Index2,Index3,Seq,reduce,lambda,w)
Output('RPseudoAAC1',pseudoaac1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRseAAC Parameters
% Calculate the type II pseudo-reduced amino acids of the protein sequences
% Index1 -> a kind of chosen physico-chemical properties of database1
% Index2 -> a kind of chosen physico-chemical properties of database1
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% lambda-> the lag of the RPseudoAAC;
% w-> the weighting factor
Index1='ARGP820103';
Index2='BEGF750102';
lambda=10;
w=0.03;
pseudoaac2=RPseudoAAC2(Index1,Index2,Seq,reduce,lambda,w);
Output('RPseudoAAC2',pseudoaac2);
%% Correlation-based features of reduced amino acids
% Parameters of the function
% Index -> a kind of chosen physico-chemical properties of database1
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% d-> the lag of the correlation 
Index='ARGP820103';
d=10;
[MBAC1,MBAC2,MBAC3]=RACF(Index,Seq,reduce,d);
Output('RACF_MBAC1',MBAC1);
Output('RACF_MBAC2',MBAC2);
Output('RACF_MBAC3',MBAC3);
%% Order-based features of reduced amino acids
% Parameters of function
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids 
% d-> the lag of the correlation
% w-> the weighting factor
d=10;
w=0.01;
[Order_d,Quasi_order] = RSeqOrder(Seq,reduce,d,w);
Output('RSeqOrder_Order_d',Order_d);
Output('RSeqOrder_Quasi_order',Quasi_order);
%% Position-based features of reduced amino acids
% Parameters of function
% Seq -> the loaded protein sequence
% reduce-> the reduced amino acids
[AAP] = RAAP(Seq,reduce)
Output('RAAP_AAP',AAP);
