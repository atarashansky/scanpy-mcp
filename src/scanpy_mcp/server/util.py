import inspect
from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
import os

from ..schema.util import *
from scmcp_shared.schema import AdataInfo
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.util import (
    forward_request,
    get_ads,
    add_op_log,
    deserialize_mcp_param,
)
from scmcp_shared.server.preset.util import ScanpyUtilMCP

logger = setup_logger()


ul_mcp = ScanpyUtilMCP().mcp
