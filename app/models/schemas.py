from __future__ import annotations
from typing import List, Dict, Optional

from pydantic import BaseModel, Field

class DERequest(BaseModel):
    datasetId: str = Field(..., example="MOCK")
    geneIds: List[str] = Field(..., example=["ENSG00000141510", "ENSG00000171862"])
    conditionA: str = Field(..., example="Primary Tumor")
    conditionB: str = Field(..., example="Solid Tissue Normal")

class DEResultRow(BaseModel):
    geneId: str
    log2FoldChange: Optional[float] = None
    pvalue: Optional[float] = None
    padj: Optional[float] = None
    baseMean: Optional[float] = None
    stat: Optional[float] = None

class DEResponse(BaseModel):
    datasetId: str
    conditionA: str
    conditionB: str
    sampleCounts: Dict[str, int]
    results: List[DEResultRow]