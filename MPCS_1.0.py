
import numpy as np
import platform

# Vector Based Calculation
def matrix_vector_cossim(matrix, vector):
    n = matrix.dot(vector)
    d = np.sqrt(matrix.sum(axis=1)) * np.sqrt(vector.sum())
    return np.max(n/d)

# Matrix Based Calculation
# switches matrix computation based on the system (can be changed to static if necessary)
def matrix_matrix_cossim(m1, m2, return_max=True):
    m1_norm = np.sqrt(m1.sum(axis=1))
    m2_norm = np.sqrt(m2.sum(axis=1))
    m = np.dot(m2, m1.T) if platform.system()[:3]=='Win' else np.matmul(m2, m1.T)
    m = m / np.outer(m2_norm, m1_norm)
    return np.max(m, axis=0), np.max(m, axis=1)

# Maximum Pairwise Cosine Similarity - MPCS
def MPCS(m1, m2, agg='mean'):
    if m1.T.shape[0] != m2.T.shape[0]:
        print('Matrices have wrong dimensions!')
        return
    
    # m1_rows * m2_rows > 2^27 elements ~ 100M matrix elements
    # matrix too large for memory => vector based calculation
    if m1.shape[0]*m2.shape[0] > 100_000_000:
        similarity1 = np.array([matrix_vector_cossim(m2, row) for row in m1])
        similarity2 = np.array([matrix_vector_cossim(m1, row) for row in m2])
    else:
        similarity1, similarity2 = matrix_matrix_cossim(m1, m2)
    
    match agg:
        case 'mean':
            agg_sim1 = np.mean(similarity1)
            agg_sim2 = np.mean(similarity2)
        case 'max':
            agg_sim1 = np.max(similarity1)
            agg_sim2 = np.max(similarity2)
        # Case 'min': Not Used
        case 'min':
            agg_sim1 = np.min(similarity1)
            agg_sim2 = np.min(similarity2)
        case _:
            print('No such aggregate function available!', agg)
            return

    symmetric_sim = max(agg_sim1, agg_sim2)
    return symmetric_sim