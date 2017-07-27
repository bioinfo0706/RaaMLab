function [m, n, S, names, sizes] = readfasta(filename)
% function [m, n, S, names, sizes] = readmfasta(filename)
% reads set of sequence from multyple fasta format.
% 'filename' is a name of a multiple fasta file;
% 'S' is an n x m matrix ('n' sequences of lenght 'm')
% 'names' is the array of strings of sequence names   
% sizes is an array of individual sequence length
% RaaMLab: a MATLAB toolbox for generating amino acid group-ings and RedAA modes
% Qi Dai, 20 Apri 2014,


MAXLENGTH = 10000000;
MAXNAME = 200;

file = fopen(filename, 'r');
disp(['Reading ',filename]);

% ----------------------------------------------------------
% now we are looking for the maximum length of the sequence

n=0;    % the number of sequences
m=0;    % the maximum length
cm = 0; % current sequence length

while 1
	[x,nr] = fscanf(file,'%c',1);
   if nr == 0 break; end;
   if x =='>'  % new sequence started
		if cm > m m=cm; end;
		cm = 0;
		fgets(file);
		n=n+1;
	else
		if isletter(x) | x=='-'| x=='('| x==')'| x=='.'
			cm=cm+1;
		end;
	end;
end

if cm > m m=cm; end;

% ----------------------------------------------------------
% go throught the file
Ss = char(m); S = [];
str = zeros(1,MAXNAME);
sizes = zeros(1,n);
frewind(file);
names=[];
i=0;j=1;
while 1
   [x,nr] = fscanf(file,'%c',1);
   if nr == 0 break; end;
	if x =='>'  % new sequence started
		if i~= 0 % save the sequence
			[x, sizes(i)]=size(Ss);
			S=strvcat(S,Ss);
			Ss = []; Ss = char(m);
		end;
		str=fgetl(file); % read the name, we remove the '>' symbol
		names=strvcat(names,str);
		disp(str);
		i=i+1;
		j=1;
	else
		if isletter(x) | x=='-'| x=='('| x==')'| x=='.'
		% processing the sequence symbol
	   	Ss(j) = upper(x); 
			j=j+1;
		end;
	end;
end

S=strvcat(S,Ss);
[x, sizes(i)]=size(Ss);

