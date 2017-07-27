function [G] = AAreduce(Index,para_dis,para_meth,maxgroup)

% Reduce amino acids into specific groups with size R according to the
% user¡¯s setting. 
% Index-> the index of databases in RaaMlab toolbox
%  para_dis-> computes the distance of amino acids based on the index of databases
%
%       'euclidean'   - Euclidean distance
%       'seuclidean'  - Standardized Euclidean distance, each coordinate
%                       in the sum of squares is inverse weighted by the
%                       sample variance of that coordinate
%       'cityblock'   - City Block distance
%       'mahalanobis' - Mahalanobis distance
%       'minkowski'   - Minkowski distance with exponent 2
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%
% para_meth-> creates a hierarchical cluster tree, using the single
%   linkage algorithm of Matlab
%      'single'    --- nearest distance
%      'complete'  --- furthest distance
%      'average'   --- unweighted average distance (UPGMA) (also known as
%                      group average)
%      'weighted'  --- weighted average distance (WPGMA)
%      'centroid'  --- unweighted center of mass distance (UPGMC) (*)
%      'median'    --- weighted center of mass distance (WPGMC) (*)
%      'ward'      --- inner squared distance (min variance algorithm) (*)
% maxgroup_> is the size of reduced amino acid set. Maxgroup can be from 1 to 20.
%
%If you want to use published amino acid groupings in the database, in which we have collected 34 kinds of published amino acid groupings. 
%You can also used this function reduce=AAreduce('index')
% Index-> the index of published amino acids in the database 4 of RaaMLab
%If you want to select no reduced method and then skip this step. You can also used this function
% reduce=AAreduce('FullAlphabet')
%   RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
%   Qi Dai, 20 Apri 2014, 

indextxt=importdata('AAindex.txt');
if nargin==1
    si=1;
    c=0;
    G=zeros(1,20);
 while c~=1 && si<=length(indextxt.data)
    a=indextxt.textdata{si};
    c=strcmp(Index,a);
    G=indextxt.data(si,:);
    si=si+1;
 end
end
if nargin==4
    for si=1:length(indextxt.textdata)
        a=indextxt.textdata{si};
        b=strcmp(Index,a);
        if b==1
            AAindex=indextxt.data(si,:);
        end
    end
    distance=pdist(AAindex',para_dis);
    Z = linkage(distance,para_meth);
    G = cluster(Z,'maxclust',maxgroup);
    G=G';
end