from fastmcp import FastMCP, Context
import os
import scanpy as sc
from ..schema.tl import *
from ..util import filter_args, add_op_log
from ..logging_config import setup_logger
logger = setup_logger(log_file=os.environ.get("SCANPYMCP_LOG_FILE", None))

tl_mcp = FastMCP("ScanpyMCP-TL-Server")


@tl_mcp.tool()
def tsne(request: TSNEModel, ctx: Context):
    """t-distributed stochastic neighborhood embedding (t-SNE) for visualization"""
    func_kwargs = filter_args(request, sc.tl.tsne)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.tsne(adata, **func_kwargs)
    add_op_log(adata, sc.tl.tsne, func_kwargs)
    return adata

@tl_mcp.tool()
def umap(request: UMAPModel, ctx: Context):
    """Uniform Manifold Approximation and Projection (UMAP) for visualization"""
    func_kwargs = filter_args(request, sc.tl.umap)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.umap(adata, **func_kwargs)
    add_op_log(adata, sc.tl.umap, func_kwargs)
    return adata

@tl_mcp.tool()
def draw_graph(request: DrawGraphModel, ctx: Context):
    """Force-directed graph drawing"""
    func_kwargs = filter_args(request, sc.tl.draw_graph)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.draw_graph(adata, **func_kwargs)
    add_op_log(adata, sc.tl.draw_graph, func_kwargs)
    return adata

@tl_mcp.tool()
def diffmap(request: DiffMapModel, ctx: Context):
    """Diffusion Maps for dimensionality reduction"""
    func_kwargs = filter_args(request, sc.tl.diffmap)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.diffmap(adata, **func_kwargs)
    adata.obsm["X_diffmap"] = adata.obsm["X_diffmap"][:,1:]
    add_op_log(adata, sc.tl.diffmap, func_kwargs)
    return adata

@tl_mcp.tool()
def embedding_density(request: EmbeddingDensityModel, ctx: Context):
    """Calculate the density of cells in an embedding"""
    func_kwargs = filter_args(request, sc.tl.embedding_density)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.embedding_density(adata, **func_kwargs)
    add_op_log(adata, sc.tl.embedding_density, func_kwargs)
    return adata

@tl_mcp.tool()
def leiden(request: LeidenModel, ctx: Context):
    """Leiden clustering algorithm for community detection"""
    func_kwargs = filter_args(request, sc.tl.leiden)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.leiden(adata, **func_kwargs)
    add_op_log(adata, sc.tl.leiden, func_kwargs)
    return adata

@tl_mcp.tool()
def louvain(request: LouvainModel, ctx: Context):
    """Louvain clustering algorithm for community detection"""
    func_kwargs = filter_args(request, sc.tl.louvain)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.louvain(adata, **func_kwargs)
    add_op_log(adata, sc.tl.louvain, func_kwargs)
    return adata

@tl_mcp.tool()
def dendrogram(request: DendrogramModel, ctx: Context):
    """Hierarchical clustering dendrogram"""
    func_kwargs = filter_args(request, sc.tl.dendrogram)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.dendrogram(adata, **func_kwargs)
    add_op_log(adata, sc.tl.dendrogram, func_kwargs)
    return adata

@tl_mcp.tool()
def dpt(request: DPTModel, ctx: Context):
    """Diffusion Pseudotime (DPT) analysis"""
    func_kwargs = filter_args(request, sc.tl.dpt)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.dpt(adata, **func_kwargs)
    add_op_log(adata, sc.tl.dpt, func_kwargs)
    return adata

@tl_mcp.tool()
def paga(request: PAGAModel, ctx: Context):
    """Partition-based graph abstraction"""
    func_kwargs = filter_args(request, sc.tl.paga)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.paga(adata, **func_kwargs)
    add_op_log(adata, sc.tl.paga, func_kwargs)
    return adata

@tl_mcp.tool()
def ingest(request: IngestModel, ctx: Context):
    """Map labels and embeddings from reference data to new data"""
    func_kwargs = filter_args(request, sc.tl.ingest)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.ingest(adata, **func_kwargs)
    add_op_log(adata, sc.tl.ingest, func_kwargs)
    return adata

@tl_mcp.tool()
def rank_genes_groups(request: RankGenesGroupsModel, ctx: Context):
    """Rank genes for characterizing groups, for differentially expressison analysis"""
    func_kwargs = filter_args(request, sc.tl.rank_genes_groups)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.rank_genes_groups(adata, **func_kwargs)
    add_op_log(adata, sc.tl.rank_genes_groups, func_kwargs)
    return adata

@tl_mcp.tool()
def filter_rank_genes_groups(request: FilterRankGenesGroupsModel, ctx: Context):
    """Filter out genes based on fold change and fraction of genes"""
    func_kwargs = filter_args(request, sc.tl.filter_rank_genes_groups)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.filter_rank_genes_groups(adata, **func_kwargs)
    add_op_log(adata, sc.tl.filter_rank_genes_groups, func_kwargs)
    return adata

@tl_mcp.tool()
def marker_gene_overlap(request: MarkerGeneOverlapModel, ctx: Context):
    """Calculate overlap between data-derived marker genes and reference markers"""
    func_kwargs = filter_args(request, sc.tl.marker_gene_overlap)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.marker_gene_overlap(adata, **func_kwargs)
    add_op_log(adata, sc.tl.marker_gene_overlap, func_kwargs)
    return adata

@tl_mcp.tool()
def score_genes(request: ScoreGenesModel, ctx: Context):
    """Score a set of genes based on their average expression"""
    func_kwargs = filter_args(request, sc.tl.score_genes)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.score_genes(adata, **func_kwargs)
    add_op_log(adata, sc.tl.score_genes, func_kwargs)
    return adata

@tl_mcp.tool()
def score_genes_cell_cycle(request: ScoreGenesCellCycleModel, ctx: Context):
    """Score cell cycle genes and assign cell cycle phases"""
    func_kwargs = filter_args(request, sc.tl.score_genes_cell_cycle)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.tl.score_genes_cell_cycle(adata, **func_kwargs)
    add_op_log(adata, sc.tl.score_genes_cell_cycle, func_kwargs)
    return adata