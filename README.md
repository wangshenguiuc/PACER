## PACER
The code is written in MATLAB. 

Please run DCA_embedding in the src folder

We need two input files to run the code NetList.

### Network file
The first one is a network file (e.g., PPI network, co-expression network). An example is the 'Network.txt' in the data folder.

Format:
>202	2173	0.588
>2173	202	0.588
>204	901	0.001
>901	204	0.001

In each network file, each row contains three numbers, representing an edge in the network. The first two numbers are gene IDs and the third is the edge weight. 

### Pathway-gene association file
The second file is a pathway-gene association file. An example is the 'Pathway_property.txt' in the data folder.

Format:
1487	2 1
3713	2 1
4201	2 1
...

The first number is the pathway ID, and the second number is the gene ID. The third number is the edge weight.

Usage: Please put the network file and pathway association file in ./data/
Running DCA_embedding.m will generate embedding vectors in ./results/ for each gene.

### Final pathway rankings.
Then, please run the Python script to get the final pathway ranking based on the embedding vector.

```bash
python embedding_top_pathways.py [-h] -f EMBED_FNAME -p PATH_FNAME -i PATH_IDX_FNAME -g GENE_IDX_FNAME
```

The tool will soon be integrated into UIUC KnowEng interface.
