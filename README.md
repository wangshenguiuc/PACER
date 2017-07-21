## PACER
The code is written in MATLAB. 

Please run DCA_embedding in the src folder

We need two input files to run the code NetList.

### Network file
The first one is a network file (e.g., PPI network, co-expression network). An example is the 'Network.txt' in the data folder.

Format:
>202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2173&nbsp;&nbsp;&nbsp;&nbsp;0.588<br />
>2173&nbsp;&nbsp;&nbsp;&nbsp;202&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.588<br />
>204&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;901&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.001<br />
>901&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;204&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0.001<br />

In each network file, each row contains three numbers, representing an edge in the network. The first two numbers are gene IDs and the third is the edge weight. 

### Pathway-gene association file
The second file is a pathway-gene association file. An example is the 'Pathway_property.txt' in the data folder.

Format:
>1487&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br />
>3713&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br />
>4201&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br />

The first number is the pathway ID, and the second number is the gene ID. The third number is the edge weight.

Usage: Please put the network file and pathway association file in ./data/
Running DCA_embedding.m will generate embedding vectors in ./results/ for each gene.

### Final pathway rankings.
Then, please run the Python script to get the final pathway ranking based on the embedding vector.

```bash
python embedding_top_pathways.py [-h] -f EMBED_FNAME -p PATH_FNAME -i PATH_IDX_FNAME -g GENE_IDX_FNAME
```

The tool will soon be integrated into UIUC KnowEng interface.
