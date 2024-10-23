import pandas
import umap
import numpy
from utilsmeta.my_logger import My_logger


logger = My_logger(log_filename = snakemake.log[0], logger_name = "process_umap")
logger = logger.get_logger()

# IDs
ids = None
with open(snakemake.input[0], 'r') as ifile:
    ids = ifile.readline().rstrip('\n').split('\t')[1:]
    logger.info('Distance matrix: {num} x {num}'.format(num=len(ids)))

# distance matrix
logger.info('Start loading distances...')
dist = numpy.loadtxt(
    fname=snakemake.input[0],
    skiprows=1,
    usecols=range(1,len(ids) + 1) # skip 1st (contains IDs)
)

# embedding
logger.info('Start embedding...')
embedding = umap.UMAP(
    n_neighbors=snakemake.params.neighbors,
    n_components=snakemake.params.components,
    min_dist=snakemake.params.min_dist,
    init='random',
    metric='precomputed',
    random_state=42
).fit_transform(dist)
logger.info('Done.')

# save to file
embedding = pandas.DataFrame(
    embedding,
    columns=['UMAP_D1', 'UMAP_D2'],
    index=ids
)
embedding.to_csv(
    path_or_buf=snakemake.output[0],
    sep='\t',
    header=True,
    index=True,
    index_label='NUCCORE_ACC'
)