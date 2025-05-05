import pytest
from fastmcp import Client
import anndata
import numpy as np
from pathlib import Path


@pytest.mark.asyncio 
async def test_subset_cells(mcp_config):
    # Pass the server directly to the Client constructor
    test_dir = Path(__file__).parent / "data/hg19"
    async with Client(mcp_config) as client:
        # First load the data
        result = await client.call_tool("io_read", {"request":{"filename": test_dir}})
        assert "AnnData" in result[0].text
        
        # Test filtering with min_genes parameter
        result = await client.call_tool("pp_subset_cells", {"request":{"min_genes": 200}})
        assert "AnnData" in result[0].text
