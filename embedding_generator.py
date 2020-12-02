import numpy as np
from minorminer import find_embedding
from dwave.system import EmbeddingComposite, DWaveSampler
import dimod

def make_complete_graph(size):
    jmatr = {}
    for i in range(size):
        for j in range(i+1, size):
            jmatr[(i,j)] = 1
    return jmatr

def count_physical_qubits(embeding):
    nq = 0
    for sub_embed in embeding:
        nq += sub_embed
    return nq

def generate_embedding(size):
    qpusampler = DWaveSampler(solver='Advantage_system1.1')
    edges = qpusampler.structure.edgelist
    gen_graph = make_complete_graph(size)
    embeddings = find_embedding(gen_graph, edges)
    n_qb = count_physical_qubits(embeddings)
    print("Embedding on %d spin variables" %(n_qb))

    return embeddings

