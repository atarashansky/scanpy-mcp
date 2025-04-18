import inspect
from fastmcp import FastMCP, Context
import os
import anndata as ad
from ..schema.util import *
from ..util import add_op_log
from ..logging_config import setup_logger
logger = setup_logger(log_file=os.environ.get("SCANPYMCP_LOG_FILE", None))


ul_mcp = FastMCP("ScanpyMCP-Util-Server")


@ul_mcp.tool()
def mark_var(request: MarkVarModel, ctx: Context):
    """
    Determine if each gene meets specific conditions and store results in adata.var as boolean values.
    For example: mitochondrion genes startswith MT-.
    The tool should be called first when calculate quality control metrics for mitochondrion, ribosomal, harhemoglobin genes, or other qc_vars.
    """
    adata = ctx.session.adata_dic[ctx.session.active_id]
    var_name = request.var_name
    gene_class = request.gene_class
    pattern_type = request.pattern_type
    patterns = request.patterns
    
    if gene_class is not None:
        if gene_class == "mitochondrion":
            adata.var["mt"] = adata.var_names.str.startswith(('MT-', 'Mt','mt-'))
            var_name = "mt"
        elif gene_class == "ribosomal":
            adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
            var_name = "ribo"
        elif gene_class == "hemoglobin":
            adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", case=False)
            var_name = "hb"
    
    if pattern_type is not None and patterns is not None:
        if pattern_type == "startswith":
            adata.var[var_name] = adata.var_names.str.startswith(patterns)
        elif pattern_type == "endswith":
            adata.var[var_name] = adata.var_names.str.endswith(patterns)
        elif pattern_type == "contains":
            adata.var[var_name] = adata.var_names.str.contains(patterns)
        else:
            raise ValueError(f"Did not support pattern_type: {pattern_type}")
    
    add_op_log(adata, mark_var, {
        "var_name": var_name,
        "gene_class": gene_class,
        "pattern_type": pattern_type,
        "patterns": patterns
    })
    
    return {var_name: adata.var[var_name].value_counts().to_dict(), "msg": f"add '{var_name}' column in adata.var"}


@ul_mcp.tool()
def list_var(request: ListVarModel, ctx: Context):
    """List key columns in adata.var. It should be called for checking when other tools need var key column names as input."""
    adata = ctx.session.adata_dic[ctx.session.active_id]
    columns = list(adata.var.columns)
    add_op_log(adata, list_var, {})
    return columns


@ul_mcp.tool()
def list_obs(request: ListObsModel, ctx: Context):
    """List key columns in adata.obs. It should be called before other tools need obs key column names input."""
    adata = ctx.session.adata_dic[ctx.session.active_id]
    columns = list(adata.obs.columns)
    add_op_log(adata, list_obs, {})
    return columns

@ul_mcp.tool()
def check_gene(request: VarNamesModel, ctx: Context):
    """Check if genes exist in adata.var_names. This tool should be called before gene expression visualizations or color by genes."""
    adata = ctx.session.adata_dic[ctx.session.active_id]
    var_names = request.var_names
    result = {v: v in adata.var_names for v in var_names}
    add_op_log(adata, check_gene, {"var_names": var_names})
    return result


@ul_mcp.tool()
def merge_adata(request: ConcatAdataModel, ctx: Context):
    """Merge multiple adata objects."""
    try:
        kwargs = {k: v for k, v in request.model_dump().items() if v is not None}
        merged_adata = ad.concat(list(ctx.session.adata_dic.values()), **kwargs)
        ctx.session.adata_dic = {}
        ctx.session.active_id = "merged_adata"
        ctx.session.adata_dic[ctx.session.active_id] = merged_adata
        
        return {"status": "success", "message": "Successfully merged all AnnData objects"}
    except Exception as e:
        logger.error(f"Error merging AnnData objects: {e}")
        raise e
