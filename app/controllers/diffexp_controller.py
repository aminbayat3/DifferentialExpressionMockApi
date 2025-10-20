from fastapi import APIRouter, HTTPException

from app.models.schemas import DERequest, DEResponse
from app.services.mock_data import generate_mock_dataset
from app.services.de_core import build_subset, run_deseq2, pack_response

router = APIRouter()

@router.post("/diffexp", response_model=DEResponse)
def differential_expression(payload: DERequest) -> DEResponse:
    dataset_id = payload.datasetId.strip()
    geneIds = [g.strip() for g in payload.geneIds if g and g.strip()]
    condition_a = payload.conditionA.strip()
    condition_b = payload.conditionB.strip()

    if not dataset_id or not geneIds or not condition_a or not condition_b:
        raise HTTPException(status_code=400, detail="Missing required fields.")

    if dataset_id.upper() == "MOCK":
        from app.services.mock_data import MOCK_GENES
        pad_n = min(1000, len(MOCK_GENES)) or 0
        padded_genes = list(set(MOCK_GENES[:pad_n]) | set(geneIds))

        counts_matrix, pheno = generate_mock_dataset(
            padded_genes, condition_a, condition_b, n_per_group=80
        )
        counts_df, clinical_df = build_subset(
            counts_matrix, pheno, padded_genes, condition_a, condition_b
        )
    else:
        raise HTTPException(status_code=400, detail="Only MOCK mode is enabled. Use datasetId='MOCK'.")

    res_df = run_deseq2(counts_df, clinical_df, ref_level="B")
    res_df = res_df.loc[res_df.index.intersection(geneIds)]

    return pack_response(
        dataset_id, condition_a, condition_b, clinical_df, res_df, keep_genes_order=geneIds
    )