import asyncio
from fastmcp import FastMCP

from .io import io_mcp
from .pp import pp_mcp
from .tl import tl_mcp
from .pl import pl_mcp
from .util import ul_mcp


scanpy_mcp = FastMCP("Scanpy-MCP-Server")


async def setup():
    await scanpy_mcp.import_server("io", io_mcp)
    await scanpy_mcp.import_server("pp", pp_mcp)
    await scanpy_mcp.import_server("tl", tl_mcp) 
    await scanpy_mcp.import_server("pl", pl_mcp) 
    await scanpy_mcp.import_server("ul", ul_mcp)