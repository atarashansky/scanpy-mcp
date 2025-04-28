
import os
import inspect
import scanpy as sc
from fastmcp import FastMCP , Context
from ..schema.pp import *
from ..util import filter_args, add_op_log
from ..logging_config import setup_logger
logger = setup_logger(log_file=os.environ.get("SCANPYMCP_LOG_FILE", None))


pp_mcp = FastMCP("ScanpyMCP-PP-Server")


@pp_mcp.tool()
def subset_cells(request: SubsetCellModel, ctx: Context):
    """subset or slice or filter cells based on total genes expressed counts and numbers. or values in adata.obs[obs_key]"""
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    func_kwargs = filter_args(request, sc.pp.filter_cells)
    if func_kwargs:
        sc.pp.filter_cells(adata, **func_kwargs)
        add_op_log(adata, sc.pp.filter_cells, func_kwargs)
    # Subset based on obs (cells) criteria
    if request.obs_key is not None:
        if request.obs_key not in adata.obs.columns:
            raise ValueError(f"Key '{request.obs_key}' not found in adata.obs")        
        mask = True  # Start with all cells selected
        if request.obs_value is not None:
            mask = mask & (adata.obs[request.obs_key] == request.obs_value)
        if request.obs_min is not None:
            mask = mask & (adata.obs[request.obs_key] >= request.obs_min)        
        if request.obs_max is not None:
            mask = mask & (adata.obs[request.obs_key] <= request.obs_max)        
        adata = adata[mask, :]
        add_op_log(adata, "subset_cells", 
            {
            "obs_key": request.obs_key, "obs_value": request.obs_value, 
             "obs_min": request.obs_min, "obs_max": request.obs_max
             }
        )
    ctx.session.adata_dic[ctx.session.active_id] = adata
    return adata


@pp_mcp.tool()
def subset_genes(request: SubsetGeneModel, ctx: Context):
    """subset/slice/filter genes based on number of cells or counts, or values in adata.var[var_key] or subset highly variable genes"""    
    func_kwargs = filter_args(request, sc.pp.filter_genes)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    if func_kwargs:
        sc.pp.filter_genes(adata, **func_kwargs)
        add_op_log(adata, sc.pp.filter_genes, func_kwargs)
    if request.var_key is not None:
        if request.var_key not in adata.var.columns:
            raise ValueError(f"Key '{request.var_key}' not found in adata.var")
        mask = True  # Start with all genes selected
        if request.var_min is not None:
            mask = mask & (adata.var[request.var_key] >= request.var_min)
        if request.var_max is not None:
            mask = mask & (adata.var[request.var_key] <= request.var_max)        
        adata = adata[:, mask]
        if request.highly_variable is not None:
            adata = adata[:, adata.var.highly_variable]
        add_op_log(adata, "subset_genes", 
            {
            "var_key": request.var_key, "var_value": request.var_value, 
             "var_min": request.var_min, "var_max": request.var_max, "hpv":  request.highly_variable
             }
        )
    ctx.session.adata_dic[ctx.session.active_id] = adata
    return adata


@pp_mcp.tool()
def calculate_qc_metrics(request: CalculateQCMetrics, ctx: Context):
    """Calculate quality control metrics(common metrics: total counts, gene number, percentage of counts in ribosomal and mitochondrial) for AnnData."""
    func_kwargs = filter_args(request, sc.pp.calculate_qc_metrics)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    func_kwargs["inplace"] = True
    sc.pp.calculate_qc_metrics(adata, **func_kwargs)
    add_op_log(adata, sc.pp.calculate_qc_metrics, func_kwargs)
    return adata


@pp_mcp.tool()
def log1p(request: Log1PModel, ctx: Context):
    """Logarithmize the data matrix (X = log(X + 1))"""
    func_kwargs = filter_args(request, sc.pp.log1p)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.log1p(adata, **func_kwargs)
        adata.raw = adata.copy()
        add_op_log(adata, sc.pp.log1p, func_kwargs)
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
        add_op_log(adata, sc.pp.normalize_total, func_kwargs)
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
    add_op_log(adata, sc.pp.pca, func_kwargs)
    return adata


@pp_mcp.tool()
def highly_variable_genes(request: HighlyVariableGenesModel, ctx: Context):
    """Annotate highly variable genes"""
    func_kwargs = filter_args(request, sc.pp.highly_variable_genes)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.pp.highly_variable_genes(adata, **func_kwargs)
    add_op_log(adata, sc.pp.highly_variable_genes, func_kwargs)
    return adata


@pp_mcp.tool()
def regress_out(request: RegressOutModel, ctx: Context):
    """Regress out (mostly) unwanted sources of variation."""
    func_kwargs = filter_args(request, sc.pp.regress_out)
    adata = ctx.session.adata_dic[ctx.session.active_id].copy()
    try:
        sc.pp.regress_out(adata, **func_kwargs)
        add_op_log(adata, sc.pp.regress_out, func_kwargs)
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
        add_op_log(adata, sc.pp.scale, func_kwargs)
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
        add_op_log(adata, sc.pp.combat, func_kwargs)
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
    add_op_log(adata, sc.pp.scrublet, func_kwargs)
    return adata


@pp_mcp.tool()
def neighbors(request: NeighborsModel, ctx: Context):
    """Compute nearest neighbors distance matrix and neighborhood graph"""
    func_kwargs = filter_args(request, sc.pp.neighbors)
    adata = ctx.session.adata_dic[ctx.session.active_id]
    sc.pp.neighbors(adata, **func_kwargs)
    add_op_log(adata, sc.pp.neighbors, func_kwargs)
    return adata
