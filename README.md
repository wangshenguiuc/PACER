# NetPath
The code is written in MATLAB. 
We need two input to run the code NetList

Format:
202	2173	0.588
2173	202	0.588
204	901	0.001
901	204	0.001
...

In each network file,  every row contains three numbers representing one edge in the network. The first two numbers are the gene id and the last number is the edge weight. 

pathway_file

Format:
1487	2
3713	2
4201	2
...

The first number is the pathway id, and the second number is the gene id.

Usage: modify the NetList and pathway_file in DCA_embedding.m
NetList is the path of network file. pathway_file is the path of pathway file.

Output will be in the same location.

The tool will be integrated into UIUC KnowEng interface soon.
