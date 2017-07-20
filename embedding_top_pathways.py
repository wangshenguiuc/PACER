### Author: Edward Huang

# from collections import OrderedDict
import file_operations
import numpy as np
import operator
import os
from scipy import linalg
import sys
import time

### For each drug, we take the top K most correlated genes, as computed by the
### drug_pathway_fisher_correlation.py script. Then, for each drug-pathway pair,
### we compute a NetPath score by multiplying each gene's correlation to that 
### drug with the gene-pathway cosine similarity scores.
### Run time: 6 minutes per embedding file.

# def create_dimension_list(network):
#     '''
#     Creates the list of dimensions for which embedding was run. PPI was run open
#     twice as many dimensions.
#     '''
#     dimension_list = [50, 100, 500]
#     if network == 'ppi':
#         # We ran more dimensions for the ppi network.
#         dimension_list += [1000, 1500, 2000]
#     return map(str, dimension_list)

def get_embedding_dct(dimension, suffix, fraction, emb_node_lst, all_genes,
    path_to_gene_dct):
    '''
    Returns a dictionary.
    Key: node (either a gene or a pathway) -> str
    Value: key's corresponding embedding vector -> list(float)
    '''
    # node_embedding_dct = OrderedDict({})
    node_embedding_dct = {}

    # First, process the input network filename.    
    # extension = '%s_0.8.%s' % (dimension, suffix)
    # if network == 'ppi':
    #     filename = '%sppi_6_net_%s' % (data_folder, extension)
    # else:
    #     filename = '%s%s.network_net_%s' % (data_folder, network, extension)

    # TODO: for PPI file.
    if isPpi:
        data_folder = './data/embedding' # With PPI.
        filename = '%s/ppi_6_net_%s_%s.%s' % (data_folder, dimension, fraction,
            suffix)
    else:
        data_folder = './embedding' # string_experimental
        filename = '%s/string_experimental_%s_%s.%s' % (data_folder, dimension,
        fraction, suffix)

    f = open(filename, 'r')
    for i, line in enumerate(f):
        # Each row in the file maps to the row in gene_pathway_id.txt.
        node = emb_node_lst[i]
        # Skip genes and pathways that aren't in expression or NCI pathways.
        if (node not in all_genes) and (node not in path_to_gene_dct):
            continue
        # Convert lines to floats.
        assert node not in node_embedding_dct
        node_embedding_dct[node] = map(float, line.split())
    f.close()
    return node_embedding_dct

# def create_pathway_vector_matrix(nci_pathways, node_embedding_dct):
#     pathway_vector_matrix = []
#     for pathway in nci_pathways:
#         pathway_vec = node_embedding_dct[pathway]
#         pathway_vector_matrix += [pathway_vec]
#     return np.matrix(pathway_vector_matrix)

# def create_gene_vector_matrix(expression_genes, node_embedding_dct):
#     gene_vector_matrix = []
#     for gene in expression_genes:
#         gene_vec = node_embedding_dct[gene]
#         gene_vector_matrix += [gene_vec]
#     return np.matrix(gene_vector_matrix)

def create_embedding_matrix(node_list, node_embedding_dct):
    '''
    Given a list of nodes (either expression genes or pathways), output a
    matrix where each row is a node's embedding vector.
    '''
    embedding_matrix = []
    for node in node_list:
        embedding_matrix += [node_embedding_dct[node]]
    return np.matrix(embedding_matrix)

def compute_drug_path_score(drug, gene_p_val_dct, node_embedding_dct,
    path_to_gene_dct):
    '''
    Given a dictionary where keys are genes and values are p-values of the
    correlation between the key and the input drug, we return a dictionary.
    Key: (drug, pathway) pair -> (str, str)
    Value: NetPath score, the sum of the cosine between each path and each gene
        of the highly correlated gene set, multiplied by that gene's p-value for
        the input drug -> float
    '''
    drug_path_score_dct = {}

    # Get embedding matrices for pathways and genes.
    nci_pathways = path_to_gene_dct.keys()
    expr_genes = gene_p_val_dct.keys()
    path_emb_matrix = create_embedding_matrix(nci_pathways, node_embedding_dct)
    gene_emb_matrix = create_embedding_matrix(expr_genes, node_embedding_dct)

    # Compute the cosine matrix between pathways and genes.
    cosine_matrix = (path_emb_matrix * gene_emb_matrix.T) / linalg.norm(
        path_emb_matrix, axis=1)[:, np.newaxis] / linalg.norm(gene_emb_matrix,
        axis=1)

    # begin TODO
    # Normalize each column of cosine matrix.
    # cosine_matrix = cosine_matrix / cosine_matrix.sum(axis=0)
    # Normalize each row of cosine matrix.
    # cosine_matrix = np.array(cosine_matrix)
    # cosine_matrix = cosine_matrix / cosine_matrix.sum(axis=1)[:,None]
    # end TODO

    for path_i, pathway in enumerate(nci_pathways):
        # Initialize the score.
        drug_path_score = 0.0
        for gene_i, gene in enumerate(expr_genes):
            # Correlation between gene and drug response.
            gene_drug_corr = gene_p_val_dct[gene]
            cos = cosine_matrix[path_i, gene_i]
            # TODO: ./debug/ppi_top_pathways_50_0.8.U_top_250_new.txt
            # drug_path_score += abs(cos * gene_drug_corr)
            drug_path_score += cos * gene_drug_corr
            # drug_path_score += abs(cos) * gene_drug_corr
            # drug_path_score += cos * abs(gene_drug_corr)
            # TODO: testing.
            # drug_path_score += cos
        drug_path_score_dct[(drug, pathway)] = drug_path_score
    return drug_path_score_dct

def write_top_pathway_file(drug_path_score_dct, results_folder, out_fname):
    # TODO: Overriding outfile. Writing in debug folder.
    # out = open(results_folder + out_fname, 'w')
    print 'writing out to debug folder...'
    out = open('./debug/test.txt', 'w')
    # Description of method.
    out.write('cosine * gene_drug_corr, unnormalized')
    out.write('\ndrug\tpath\tscore\n')
    for (drug, pathway), score in drug_path_score_dct:
        out.write('%s\t%s\t%f\n' % (drug, pathway, score))
    out.close()

# Writes out to file the top drug-pathway pairs with score equal to the inverse
# of their overall ranking. 
def write_inverse_rankings(results_folder, filename):   
    brd_drug_to_name_dct = file_operations.get_brd_drug_to_name_dct()

    drug_path_dct = {}
    f = open(results_folder + filename, 'r')
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        line = line.strip().split('\t')
        assert len(line) == 3
        drug, path, score = line
        drug_path_dct[(drug, path)] = float(score)
    f.close()

    # Sort ppi dictionary by value.
    ranked_drug_path_dct = {}
    drug_path_dct = sorted(drug_path_dct.items(), key=operator.itemgetter(1),
        reverse=True)
    for i, ((drug, path), score) in enumerate(drug_path_dct):
        inverse_ranking = (i + 1) / float(len(drug_path_dct))
        ranked_drug_path_dct[(drug, path)] = inverse_ranking

    ranked_drug_path_dct = sorted(ranked_drug_path_dct.items(),
        key=operator.itemgetter(1), reverse=True)

    inverse_folder = '%sinverse_rankings/' % results_folder
    if not os.path.exists(inverse_folder):
        os.makedirs(inverse_folder)

    out = open(inverse_folder + 'inverse_' + filename, 'w')
    out.write('drug\tpath\tinverse_rank\n')
    for (drug, path), score in ranked_drug_path_dct:
        out.write('%s\t%s\t%g\n' % (brd_drug_to_name_dct[drug], path, score))
    out.close()

# def compute_drug_pathway_scores(network, top_k):
def compute_drug_pathway_scores():
    # Extract the NCI pathway data.
    path_to_gene_dct, nci_genes = file_operations.get_path_to_gene_dct()
    # Find genes and pathways that appear in embedding.
    emb_node_lst = file_operations.get_emb_node_lst()
    # Find the most significantly correlated genes for each drug.
    all_genes, drug_corr_genes_dct = file_operations.get_drug_corr_genes_dct(
        top_k, emb_node_lst)

    # TODO
    dimension_list = map(str, [50, 100, 500, 1000])
    fraction_list, suffix_list = map(str, [0.3, 0.5, 0.8]), ['U', 'US']
    if isPpi:
        dimension_list, fraction_list, suffix_list = ['50'], ['0.8'], ['U']

    for dimension in dimension_list:
        for fraction in fraction_list:
            for suffix in suffix_list:
                # Get the embedding vectors for each node (gene or pathway).
                node_embedding_dct = get_embedding_dct(dimension, suffix,
                    fraction, emb_node_lst, all_genes, path_to_gene_dct)

                # Calculate the score for each drug-pathway pair.
                drug_path_score_dct = {}
                for drug in drug_corr_genes_dct:
                    gene_p_val_dct = drug_corr_genes_dct[drug]
                    drug_path_score_dct.update(compute_drug_path_score(drug,
                        gene_p_val_dct, node_embedding_dct, path_to_gene_dct))

                # Sort the score dictionary, and write to file.
                drug_path_score_dct = sorted(drug_path_score_dct.items(),
                    key=operator.itemgetter(1), reverse=True)
                
                results_folder = './results/embedding/'
                if not os.path.exists(results_folder):
                    os.makedirs(results_folder)

                # out_fname = '%s_top_pathways_%s_embed%d.txt' % (network, extension,
                    # top_k)
                # TODO: for PPI file.
                if isPpi:
                    out_fname = 'top_pathways_ppi_%s_%s_%s_embed_%d.txt' % (
                        dimension, fraction, suffix, top_k)
                else:
                    out_fname = 'top_pathways_%s_%s_%s_embed_%d.txt' % (
                        dimension, fraction, suffix, top_k)
                write_top_pathway_file(drug_path_score_dct, results_folder,
                    out_fname)
                # TODO: Currently not writing out inverse rankings.
                # write_inverse_rankings(results_folder, out_fname)

def main():
    # if (len(sys.argv) != 3):
    if (len(sys.argv) not in [2, 3]):
        # print "Usage: " + sys.argv[0] + " dca_network top_k"
        print "Usage: " + sys.argv[0] + " top_k ppi<optional>"
        exit(1)
    # network = sys.argv[1]
    # assert network in ['ppi', 'genetic', 'literome', 'sequence']
    # top_k = int(sys.argv[2])
    global isPpi, top_k
    isPpi = False
    top_k = int(sys.argv[1])
    if len(sys.argv) == 3:
        isPpi = True
        assert sys.argv[2] == 'ppi'
    
    # compute_drug_pathway_scores(network, top_k)
    compute_drug_pathway_scores()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))