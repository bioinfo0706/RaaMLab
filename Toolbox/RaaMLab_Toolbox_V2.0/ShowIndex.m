function ShowIndex(database,index)
% show the databases in the RaaMlab toolbox
% database -> kinds of the databases in the RaaMlab Choices are:
%       '1'  - AAindex1 for amino acid index of 20 numerical values
%       '2'  - AAindex2 for the amino acid mutation matrix
%       '3'  - AAindex3 for the statistical protein contact potentials
%       '4'  - Published amino acid groupings
% index    -> the index of databases in RaaMlab toolbox
% if you don't know the index of databases, you can use the
% ShowIndex(database) to browe all the index in this database
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014, 
Names=caseread('Index.txt');
[s1 s2]=size(Names);
Keyindex=zeros(1,1000);
Keyindex(1)=1;k=1;
Worlds=cell(1,1000);
% the index of all the data, 1-551 for database1, 552-645 for database2,
% 646-692 for databses3, and 693-726 for databse4
KeyWorld=Names(1,3:s2);%
for i=1:s1
    if findstr('//',Names(i,:))==1
        k=k+1;
        Keyindex(k)=i;%position
        if i<=s1-1
            KeyWorld=strcat(KeyWorld,Names(i+1,3:s2));
            Worlds{k-1}=Names(i+1,:);
        end
    end
end

if nargin ==0
    error('You should input the righ parameters database and AAindex, If you don not know AAindex, you can input database only, you will find all the AAindex in this database.');
elseif nargin ==1
    if database==1|database==2|database==3|database==4
    switch(database)
    case 1
        for t=1:551
            in=['No.',num2str(t),' AAindex in the database ',num2str(database), ' is  ',Worlds{t}];
            disp(in)
        end
     case 2
        for t=552:645
            in=['No.',num2str(t-551),' AAindex of the database ',num2str(database), ' is  ',Worlds{t}];
            disp(in)
        end
     case 3
        for t=646:692
            in=['No.',num2str(t-645),' AAindex of the database ',num2str(database), ' is  ',Worlds{t}];
            disp(in)
        end
     case 4
        for t=693:726
            in=['No.',num2str(t-692),' AAindex of the database ',num2str(database), ' is  ',Worlds{t}];
            disp(in)
        end
    end
    else
      error('You should input the righ parameters database. If you don not know AAindex, you can input database only, you will find all the AAindex in this database'); 
    end
else
 dw1=0;dw=0;ii=1;
while dw1==0
   dw1=findstr(index,Worlds{ii}) ;
   if dw1~=0
   dw=ii;
   else
       dw1=0;
       ii=ii+1;
   end
end
 dw=dw+1;
   if dw~=0 
   display(Names(Keyindex(dw):Keyindex(dw+1),:));
   else
    warnning=strcat('There is no: ',index);
    display(warnning);
   end     
end

