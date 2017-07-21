import argparse
import numpy as np
import operator
from scipy import linalg

### For each drug, we get the most correlated genes. Then, for each drug-pathway
### pair, we compute a PACER score by multiplying each gene's correlation to
### that drug with the gene-pathway cosine similarity scores.

def get_embedding_matrix():
    '''
    Returns the embedding matrix, where each row is a node in the network.
    Number of columns is equal to the number of dimensions.
    '''
    embedding_matrix = []
    f = open(args.embed_fname, 'r')
    for line in f:
        embedding_matrix += [map(float, line.split())]
    f.close()
    return np.matrix(embedding_matrix)

def compute_drug_path_score(drug, gene_corr_dct, embedding_matrix, path_to_gene_dct):
    '''
    Given a dictionary where keys are genes and values are p-values of the
    correlation between the key and the input drug, we return a dictionary.
    Key: (drug, pathway) pair -> (str, str)
    Value: PACER score, the sum of the cosine between each path and each gene
        of the highly correlated gene set, multiplied by that gene's p-value for
        the input drug -> float
    '''
    # Get the pathway and gene index lists.
    pathway_lst = []
    f = open(args.path_idx_fname, 'r')
    for line in f:
        pathway_lst += [int(line.strip())]
    f.close()
    gene_lst = []
    f = open(args.gene_idx_fname, 'r')
    for line in f:
        gene_lst += [int(line.strip())]
    f.close()

    # Assert gene and pathway indices do not overlap.
    assert len(set(pathway_lst).intersection(gene_lst)) == 0
    # Fetch the embedding matrices for pathways and for genes.
    pathway_embedding_matrix = embedding_matrix[pathway_lst]
    gene_embedding_matrix = embedding_matrix[gene_lst]

    # Compute the cosine matrix between all pairs of pathways and genes.
    cosine_matrix = (pathway_embedding_matrix * gene_embedding_matrix.T) / linalg.norm(
        pathway_embedding_matrix, axis=1)[:, np.newaxis] / linalg.norm(gene_embedding_matrix,
        axis=1)

    drug_path_score_dct = {}
    for pathway in path_to_gene_dct:
        # Initialize the score.
        drug_path_score = 0.0
        for gene in gene_corr_dct:
            # Get the correlation between gene and drug response.
            gene_drug_corr = gene_corr_dct[gene]
            # Get the cosine similarity between the pathway and gene.
            cos = cosine_matrix[pathway, gene]
            # Add to the pathway the product between the cosine similarity and the correlation.
            drug_path_score += cos * gene_drug_corr
        drug_path_score_dct[(drug, pathway)] = drug_path_score
    return drug_path_score_dct

def write_top_pathway_file(drug_path_score_dct, results_folder):
    out_fname = './results/top_pathways_%s.txt' % args.embed_fname
    out = open(out_fname, 'w')
    out.write('Drug\tPathway\tscore\n')
    for (drug, pathway), score in drug_path_score_dct:
        out.write('%s\t%s\t%f\n' % (drug, pathway, score))
    out.close()

def compute_drug_pathway_scores():
    '''
    Compute the score between each drug and each pathway. Write out the scores
    for each pair.
    '''
    # Read the ordered pathway dictionary mapping pathways to gene sets.
    path_to_gene_dct = {}
    f = open(args.path_fname, 'r')
    for line in f:
        pathway, gene, weight = line.strip()
        if pathway not in path_to_gene_dct:
            path_to_gene_dct[pathway] = set([])
        path_to_gene_dct[pathway].add(gene)
    f.close()

    # Read the ordered dictionary mapping drugs to the top k correlated genes.
    # Each value is another dictionary mapping each gene to the Pearson correlation.
    with open(args.drug_fname, 'r') as fp:
        drug_to_gene_dct = json.load(fp)
    fp.close()

    # Read the embedding vectors for each node (gene or pathway).
    embedding_matrix = get_embedding_matrix()

    # Calculate the score for each drug-pathway pair.
    drug_path_score_dct = {}
    for drug in drug_to_gene_dct:
        # Get the correlation gene dictionary for the current drug.
        gene_corr_dct = drug_to_gene_dct[drug]
        # Compute PACER scores for the current drug.
        drug_path_score_dct.update(compute_drug_path_score(drug, gene_corr_dct,
            embedding_matrix, path_to_gene_dct))

    # Sort the score dictionary, and write to file.
    drug_path_score_dct = sorted(drug_path_score_dct.items(), key=operator.itemgetter(1),
        reverse=True)
    
    write_top_pathway_file(drug_path_score_dct, results_folder)

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--embed_fname', help='Input file name for embedding.',
        required=True, type=str)
    parser.add_argument('-p', '--path_fname', help='File name of json dictionary mapping pathway names to gene sets.',
        required=True, type=str)
    parser.add_argument('-d', '--drug_fname', help='File name of json dictionary mapping drugs to gene correlation dictionaries.',
        required=True, type=str)
    parser.add_argument('-i', '--path_idx_fname', help='File name of pathway indices corresponding to rows in embedding file.',
        required=True, type=str)
    parser.add_argument('-g', '--gene_idx_fname', help='File name of gene indices corresponding to rows in embedding file.',
        required=True, type=str)
    args = parser.parse_args()

def main():
    parse_args()
    compute_drug_pathway_scores()

if __name__ == '__main__':
    main()