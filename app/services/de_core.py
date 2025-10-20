from __future__ import annotations
from typing import List, Optional, Tuple

import pandas as pd
from fastapi import HTTPException

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference


def build_subset(
    counts: pd.DataFrame,   # genes x samples
    pheno: pd.DataFrame,    # rows = samples
    geneIds: List[str],
    condition_a: str,
    condition_b: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    common_samples = counts.columns.intersection(pheno.index)
    if len(common_samples) < 4:
        raise HTTPException(status_code=400, detail="Not enough overlapping samples between counts and phenotype.")

    ph = pheno.loc[common_samples].copy()
    ph = ph.loc[ph["sample_type"].isin([condition_a, condition_b])]
    if ph.empty:
        raise HTTPException(status_code=400, detail="No samples match the requested conditions.")

    cnt = counts.loc[:, ph.index].copy()
    genes_present = [g for g in geneIds if g in cnt.index]
    if len(genes_present) == 0:
        raise HTTPException(status_code=400, detail="None of the requested genes were found in the counts matrix.")

    cnt = cnt.loc[genes_present]

    counts_for_deseq = cnt.T  # samples x genes

    clinical = pd.DataFrame(
        {"condition": ph["sample_type"].replace({condition_a: "A", condition_b: "B"})},
        index=ph.index,
    )
    clinical["condition"] = pd.Categorical(clinical["condition"], categories=["B", "A"], ordered=True)

    nA = int((clinical["condition"] == "A").sum())
    nB = int((clinical["condition"] == "B").sum())
    if nA < 2 or nB < 2:
        raise HTTPException(status_code=400, detail=f"Too few samples per group (A={nA}, B={nB}).")

    return counts_for_deseq, clinical


def _filter_informative_genes_groupwise(
    counts_df: pd.DataFrame,   # samples x genes
    clinical_df: pd.DataFrame,
    min_nonzero_per_group: int,
    min_total_counts: int,
) -> pd.DataFrame:
    group = clinical_df["condition"]
    a_idx = group == "A"
    b_idx = group == "B"

    nzA = (counts_df.loc[a_idx] > 0).sum(axis=0)
    nzB = (counts_df.loc[b_idx] > 0).sum(axis=0)
    var_all = counts_df.var(axis=0)
    total_all = counts_df.sum(axis=0)

    keep = (
        (nzA >= min_nonzero_per_group) &
        (nzB >= min_nonzero_per_group) &
        (var_all > 0) &
        (total_all >= min_total_counts)
    )
    return counts_df.loc[:, keep]


def run_deseq2(
    counts_df: pd.DataFrame,  # samples x genes
    clinical_df: pd.DataFrame,
    ref_level: str = "B",
) -> pd.DataFrame:
    """Run PyDESeq2 using the documented API (design='~condition'), with robust filtering and retry."""
    counts_df = counts_df.round().astype(int)

    counts_df1 = _filter_informative_genes_groupwise(
        counts_df, clinical_df,
        min_nonzero_per_group=2, min_total_counts=10
    )
    if counts_df1.shape[1] < 10:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Too few informative genes after filtering ({counts_df1.shape[1]}). "
                f"Add more mock genes or widen your selection."
            )
        )

    inference = DefaultInference(n_cpus=1)

    try:
        dds = DeseqDataSet(
            counts=counts_df1,
            metadata=clinical_df,
            design="~condition",
            refit_cooks=True,
            inference=inference,
        )
        dds.deseq2()
        used_counts = counts_df1
    except Exception:
        counts_df2 = _filter_informative_genes_groupwise(
            counts_df1, clinical_df,
            min_nonzero_per_group=5, min_total_counts=50
        )
        if counts_df2.shape[1] < 10:
            raise HTTPException(
                status_code=400,
                detail=(
                    "DE failed even after stricter filtering. "
                    "Add more mock genes or upgrade to pydeseq2>=0.5.1."
                ),
            )
        dds = DeseqDataSet(
            counts=counts_df2,
            metadata=clinical_df,
            design="~condition",
            refit_cooks=True,
            inference=inference,
        )
        dds.deseq2()
        used_counts = counts_df2

    stat_res = DeseqStats(
        dds,
        contrast=["condition", "A", ref_level],
        inference=inference,
    )
    stat_res.summary()

    res = stat_res.results_df.copy()
    res.index.name = "geneId"
    res = res.loc[res.index.intersection(used_counts.columns)]
    return res


def pack_response(
    dataset_id: str,
    condition_a: str,
    condition_b: str,
    clinical_df: pd.DataFrame,
    res_df: pd.DataFrame,
    keep_genes_order: Optional[List[str]] = None,
):
    from app.models.schemas import DEResponse, DEResultRow  # avoid circular import
    import pandas as pd

    counts = clinical_df["condition"].replace({"A": condition_a, "B": condition_b})
    sample_counts = counts.value_counts().to_dict()

    if keep_genes_order:
        res_df = res_df.reindex(keep_genes_order)

    results = [
        DEResultRow(
            geneId=str(g),
            log2FoldChange=(row.get("log2FoldChange") if pd.notna(row.get("log2FoldChange")) else None),
            pvalue=(row.get("pvalue") if pd.notna(row.get("pvalue")) else None),
            padj=(row.get("padj") if pd.notna(row.get("padj")) else None),
            baseMean=(row.get("baseMean") if pd.notna(row.get("baseMean")) else None),
            stat=(row.get("stat") if pd.notna(row.get("stat")) else None),
        )
        for g, row in res_df.to_dict(orient="index").items()
    ]

    return DEResponse(
        datasetId=dataset_id,
        conditionA=condition_a,
        conditionB=condition_b,
        sampleCounts={k: int(v) for k, v in sample_counts.items()},
        results=results,
    )