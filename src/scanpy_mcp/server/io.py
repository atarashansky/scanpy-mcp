import os
import inspect
from pathlib import Path
import scanpy as sc
from fastmcp import FastMCP , Context
from ..schema.io import *
from ..util import filter_args
from ..logging_config import setup_logger
logger = setup_logger(log_file=os.environ.get("SCANPYMCP_LOG_FILE", None))

io_mcp = FastMCP("ScanpyMCP-IO-Server")


@io_mcp.tool()
def read(request: WriteModel, ctx: Context):
    """
    Read data from various file formats (h5ad, 10x, text files, etc.) or directory path.
    """
    kwargs = request.model_dump()
    file = Path(kwargs.get("filename", None))
    if hasattr(ctx.session, "adata_dic"):
        if kwargs.get("sampleid", None) is not None:
            ctx.session.active_id = kwargs["sampleid"]
        else:
            ctx.session.active_id = f"adata{len(ctx.session.adata_dic)}"
    else:
        ctx.session.adata_dic = {}
        ctx.session.active_id = "adata0"

    if file.is_dir():
        kwargs["path"] = kwargs["filename"]
        func_kwargs = filter_args(request, sc.read_10x_mtx)
        adata = sc.read_10x_mtx(kwargs["path"], **func_kwargs)
    elif file.is_file():
        func_kwargs = filter_args(request, sc.read)
        adata = sc.read(**func_kwargs)
        if not kwargs.get("first_column_obs", True):
            adata = adata.T
    else:
        raise ValueError("filename must be a file or a directory")
    adata.layers["counts"] = adata.X
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    ctx.session.adata_dic[ctx.session.active_id] = adata
    return adata


@io_mcp.tool()
def write(request: WriteModel, ctx: Context):
    """save adata into a file.
    """    
    kwargs = request.model_dump()
    args = request.model_fields_set
    parameters = inspect.signature(sc.write).parameters
    func_kwargs = {k: kwargs.get(k) for k in args if k in parameters}
    sc.write(func_kwargs["filename"], ctx.session.adata_dic[ctx.session.active_id])
    return {"filename": func_kwargs["filename"], "msg": "success to save file"}
