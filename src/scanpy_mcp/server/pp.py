
import os
import inspect
import scanpy as sc
from fastmcp import FastMCP , Context
from ..schema.pp import *
from ..util import filter_args
from ..logging_config import setup_logger
logger = setup_logger(log_file=os.environ.get("SCANPYMCP_LOG_FILE", None))


pp_mcp = FastMCP("ScanpyMCP-PP-Server")


@pp_mcp.tool()
def filter_cells(request: FilterCells, ctx: Context):
    """Filter cells based on counts and numbers of genes expressed."""
    func_kwargs = filter_args(request, sc.pp.filter_cells)
    logger.info(ctx.session.adata_dic[ctx.session.active_id])
    adata = ctx.session.adata_dic[ctx.session.active_id]    
    sc.pp.filter_cells(adata, **func_kwargs)
    return adata


@pp_mcp.tool()
def filter_genes(request: FilterGenes, ctx: Context):
    """Filter genes based on number of cells or counts"""    
    func_kwargs = filter_args(request, sc.pp.filter_genes)
    adata = ctx.session.adata_dic[ctx.session.active_id]    
    sc.pp.filter_genes(adata, **func_kwargs)
    return adata


@pp_mcp.tool()
def calculate_qc_metrics(request: CalculateQCMetrics, ctx: Context):
    """Calculate quality control metrics(common metrics: total counts, gene number, percentage of counts in ribosomal and mitochondrial) for AnnData."""
    func_kwargs = filter_args(request, sc.pp.calculate_qc_metrics)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    func_kwargs["inplace"] = True
    sc.pp.calculate_qc_metrics(adata, **func_kwargs)
    return adata


@pp_mcp.tool()
def log1p(request: Log1PModel, ctx: Context):
    """Logarithmize the data matrix (X = log(X + 1))"""
    func_kwargs = filter_args(request, sc.pp.log1p)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.log1p(adata, **func_kwargs)
        adata.raw = adata.copy()
    except Exception as e:
        raise e
    ctx.session.adata_dic[ctx.session.active_id] = adata
    return adata


@pp_mcp.tool()
def normalize_total(request: NormalizeTotalModel, ctx: Context):
    """Normalize counts per cell to the same total count"""
    func_kwargs = filter_args(request, sc.pp.normalize_total)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.normalize_total(adata, **func_kwargs)
    except Exception as e:
        raise e
    ctx.session.adata_dic[ctx.session.active_id] = adata
    return adata


@pp_mcp.tool()
def pca(request: PCAModel, ctx: Context):
    """Principal component analysis"""
    func_kwargs = filter_args(request, sc.pp.pca)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.pp.pca(adata, **func_kwargs)
    return adata


@pp_mcp.tool()
def highly_variable_genes(request: HighlyVariableGenesModel, ctx: Context):
    """Annotate highly variable genes"""
    func_kwargs = filter_args(request, sc.pp.highly_variable_genes)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.pp.highly_variable_genes(adata, **func_kwargs)
    return adata


@pp_mcp.tool()
def regress_out(request: RegressOutModel, ctx: Context):
    """Regress out (mostly) unwanted sources of variation."""
    func_kwargs = filter_args(request, sc.pp.regress_out)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.regress_out(adata, **func_kwargs)
    except Exception as e:
        raise e
    ctx.session.adata_dic[ctx.session.active_id] = adata
    return adata


@pp_mcp.tool()
def scale(request: ScaleModel, ctx: Context):
    """Scale data to unit variance and zero mean"""
    func_kwargs = filter_args(request, sc.pp.scale)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.scale(adata, **func_kwargs)
    except Exception as e:
        raise e
    ctx.session.adata_dic[ctx.session.active_id] = adata        
    return adata


@pp_mcp.tool()
def combat(request: CombatModel, ctx: Context):
    """ComBat function for batch effect correction"""
    func_kwargs = filter_args(request, sc.pp.combat)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.combat(adata, **func_kwargs)
    except Exception as e:
        raise e
    ctx.session.adata_dic[ctx.session.active_id] = adata         
    return adata


@pp_mcp.tool()
def scrublet(request: ScrubletModel, ctx: Context):
    """Predict doublets using Scrublet"""
    func_kwargs = filter_args(request, sc.pp.scrublet)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.pp.scrublet(adata, **func_kwargs)
    return adata


@pp_mcp.tool()
def neighbors(request: NeighborsModel, ctx: Context):
    """Compute nearest neighbors distance matrix and neighborhood graph"""
    func_kwargs = filter_args(request, sc.pp.neighbors)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.pp.neighbors(adata, **func_kwargs)
    return adata
