import os
import inspect
from functools import partial
import scanpy as sc
from fastmcp import FastMCP, Context
from ..schema.pl import *
from pathlib import Path
from ..logging_config import setup_logger
from ..util import filter_args, set_fig_path
from ..logging_config import setup_logger
logger = setup_logger(log_file=os.environ.get("SCANPYMCP_LOG_FILE", None))

pl_mcp = FastMCP("ScanpyMCP-PL-Server")



@pl_mcp.tool()
def pca(request: PCAModel, ctx: Context):
    """Scatter plot in PCA coordinates. default figure for PCA plot"""
    func_kwargs = filter_args(request, sc.pl.pca)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"

    fig = sc.pl.pca(adata, **func_kwargs)
    fig_path = set_fig_path("pca", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def diffmap(request: DiffusionMapModel, ctx: Context):
    """Plot diffusion map embedding of cells."""
    func_kwargs = filter_args(request, sc.pl.diffmap)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.diffmap(adata, **func_kwargs)
    fig_path = set_fig_path("diffmap", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def violin(request: ViolinModel, ctx: Context):
    """Plot violin plot of one or more variables."""
    func_kwargs = filter_args(request, sc.pl.violin)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.violin(adata, **func_kwargs)
    fig_path = set_fig_path("violin", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def stacked_violin(request: StackedViolinModel, ctx: Context):
    """Plot stacked violin plots. Makes a compact image composed of individual violin plots stacked on top of each other."""
    func_kwargs = filter_args(request, sc.pl.stacked_violin)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.stacked_violin(adata, **func_kwargs)
    fig_path = set_fig_path("stacked_violin", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def heatmap(request: HeatmapModel, ctx: Context):
    """Heatmap of the expression values of genes."""
    func_kwargs = filter_args(request, sc.pl.heatmap)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.heatmap(adata, **func_kwargs)
    fig_path = set_fig_path("heatmap", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def dotplot(request: DotplotModel, ctx: Context):
    """Plot dot plot of expression values per gene for each group."""
    func_kwargs = filter_args(request, sc.pl.dotplot)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.dotplot(adata, **func_kwargs)
    fig_path = set_fig_path("dotplot", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def matrixplot(request: MatrixplotModel, ctx: Context):
    """matrixplot, Create a heatmap of the mean expression values per group of each var_names."""
    func_kwargs = filter_args(request, sc.pl.matrixplot)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.matrixplot(adata, **func_kwargs)
    fig_path = set_fig_path("matrixplot", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def tracksplot(request: TracksplotModel, ctx: Context):
    """tracksplot, compact plot of expression of a list of genes."""
    func_kwargs = filter_args(request, sc.pl.tracksplot)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.tracksplot(adata, **func_kwargs)
    fig_path = set_fig_path("tracksplot", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def scatter(request: EnhancedScatterModel, ctx: Context):
    """Plot a scatter plot of two variables, Scatter plot along observations or variables axes."""
    func_kwargs = filter_args(request, sc.pl.scatter)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.scatter(adata, **func_kwargs)
    fig_path = set_fig_path("scatter", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def embedding(request: EmbeddingModel, ctx: Context):
    """Scatter plot for user specified embedding basis (e.g. umap, tsne, etc)."""
    func_kwargs = filter_args(request, sc.pl.embedding)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.embedding(adata, **func_kwargs)
    fig_path = set_fig_path("embedding", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def embedding_density(request: EmbeddingDensityModel, ctx: Context):
    """Plot the density of cells in an embedding."""
    func_kwargs = filter_args(request, sc.pl.embedding_density)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.embedding_density(adata, **func_kwargs)
    fig_path = set_fig_path("embedding_density", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def rank_genes_groups(request: RankGenesGroupsModel, ctx: Context):
    """Plot ranking of genes based on differential expression."""
    func_kwargs = filter_args(request, sc.pl.rank_genes_groups)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.rank_genes_groups(adata, **func_kwargs)
    fig_path = set_fig_path("rank_genes_groups", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def rank_genes_groups_dotplot(request: RankGenesGroupsDotplotModel, ctx: Context):
    """Plot ranking of genes(DEGs) using dotplot visualization. Defualt plot DEGs for rank_genes_groups tool"""
    func_kwargs = filter_args(request, sc.pl.rank_genes_groups_dotplot)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.rank_genes_groups_dotplot(adata, **func_kwargs)
    fig_path = set_fig_path("rank_genes_groups_dotplot", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def clustermap(request: ClusterMapModel, ctx: Context):
    """Plot hierarchical clustering of cells and genes."""
    func_kwargs = filter_args(request, sc.pl.clustermap)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.clustermap(adata, **func_kwargs)
    fig_path = set_fig_path("clustermap", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def highly_variable_genes(request: HighlyVariableGenesModel, ctx: Context):
    """plot highly variable genes; Plot dispersions or normalized variance versus means for genes."""
    func_kwargs = filter_args(request, sc.pl.highly_variable_genes)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.highly_variable_genes(adata, **func_kwargs)
    fig_path = set_fig_path("highly_variable_genes", **func_kwargs)
    return {"figpath": fig_path}


@pl_mcp.tool()
def pca_variance_ratio(request: PCAVarianceRatioModel, ctx: Context):
    """Plot the PCA variance ratio to visualize explained variance."""
    func_kwargs = filter_args(request, sc.pl.pca_variance_ratio)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    
    func_kwargs.pop("return_fig", True)
    func_kwargs["show"] = False
    func_kwargs["save"] = ".png"
    
    fig = sc.pl.pca_variance_ratio(adata, **func_kwargs)
    fig_path = set_fig_path("pca_variance_ratio", **func_kwargs)
    return {"figpath": fig_path}

