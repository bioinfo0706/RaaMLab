   RaaMLab is a free and open-source MATLAB toolbox that is to generate amino acid groupings using amino acids¡¯ properties and classification methods, and further extract the structural and physicochemical features of reduced amino acids. It offers four kinds of databases on physico-chemical properties of amino acids and amino acid groupings, 49 amino acid classification methods, and five kinds of biophyscio-chemical features of reduced amino acids including content-based features, correlation-based features, order-based features, position-based features and pseudo-reduced amino acid compositions, which can be easily computed based on the amino acid groupings with the user-defined alphabet size and amino acids¡¯ properties. 



Install RaaMLab

   RaaMLab has been successfully tested on Linux and Windows systems. The author could download the RaaMLab from http://bioinfo.zstu.edu.cn/RaaMLab.htm. The RAAM Toolbox requires at least MATLAB Release 7. No compilation is required and the use of MATLAB as a basis allows the toolbox to be used on Windows, Linux, Unix and MAC OS machines. The install process of RaaMLab is very easy:

   On Windows:
   (1): download the RaaMLab (.zip)
   (2): extract or uncompress the .zip file 
   (3): run startup.m

  On Linux:
  (1): download the propy package (.tar.gz) 
  (2): tar -zxf propy-1.0.tar.gz 
  (3): ./matlab
  (4): run startup.m 


Main functions


   startup.m   % Setup the RaaMLab toolbox
   readfasta.m % Reading testseq.txt
   ShowIndex.m % Show all the official information of the databases
   AAreduce.m  % Reduced amino acids using the methods in RaaMLab
   RAAC.m      % Calculate the k-mer content distribution of the reduced amino acids
   RCTD.m      % Compute the composition, transition and distribution of the reduced amino acids
   RPseudoAAC1.m % Calculate the type I pseudo-reduced amino acids of the protein sequences
   RPseudoAAC2.m % Calculate the type II pseudo-reduced amino acids of the protein sequences
   RACF.m        % Compute autocorrelation features, normalized Moreau¨CBroto autocorrelation, Moran autocorrelation and Geary autocorrelation
   RSeqOrder.m   % Compute two kinds of order-based features, one is sequence-order-coupling number, and the other is quasi-sequence-order.
   RAAP.m        % Compute the position-based features of reduced amino acids.


Supplementary material

   UserGuide.pdf   % A userguide for using RaaMLab toolbox including some examples for each functions
   testseq.txt     % An example of the datasets
   Index.txt       % The database file
   Grantham-distance.txt % Grantham chemical distance matrix of amino acids
   Schneider-Wrede-distance.txt % Schneider¨CWrede physicochemical distance matrix
   



If you find mistakes, or have suggestions for improvements, please send them to the mailing list: daiailiu04@yahoo.com



