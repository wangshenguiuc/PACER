PACER
The code is written in MATLAB. 

Please run DCA_embedding in the src folder

We need two input to run the code NetList. 
The first one is a network file (e.g., PPI network, co-expression network). An example is the 'Network.txt' in the data folder.

Format:
202	2173	0.588
2173	202	0.588
204	901	0.001
901	204	0.001
...

In each network file,  every row contains three numbers representing one edge in the network. The first two numbers are the gene id and the last number is the edge weight. 

The second file is a pathway gene association file. An example is the 'Pathway_property.txt' in the data folder.

Format:
1487	2 1
3713	2 1
4201	2 1
...

The first number is the pathway id, and the second number is the gene id. The third number is the weight.

Usage: Please put the network file and pathway association file in the data folder.
Run DCA_embedding will generate embedding vector for each gene in result folder.

Then please run the Python script to get the final pathway ranking based on the embedding vector.

The tool will be integrated into UIUC KnowEng interface soon.
