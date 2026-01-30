"""数据分析路由"""
from fastapi import APIRouter
from pydantic import BaseModel
from typing import Optional
from services.scfocus_service import scfocus_service

router = APIRouter()


class PreprocessRequest(BaseModel):
    min_genes: int = 200
    min_cells: int = 3
    n_top_genes: int = 2000
    n_pcs: int = 50
    n_neighbors: int = 15
    use_existing_embedding: Optional[str] = None


class PreprocessResponse(BaseModel):
    success: bool
    message: str
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    n_hvg: Optional[int] = None
    n_pcs: Optional[int] = None
    has_umap: Optional[bool] = None


@router.post("/preprocess/{session_id}", response_model=PreprocessResponse)
async def preprocess_data(session_id: str, request: PreprocessRequest):
    """预处理单细胞数据"""
    result = scfocus_service.preprocess_data(
        session_id=session_id,
        min_genes=request.min_genes,
        min_cells=request.min_cells,
        n_top_genes=request.n_top_genes,
        n_pcs=request.n_pcs,
        n_neighbors=request.n_neighbors,
        use_existing_embedding=request.use_existing_embedding
    )
    return result


@router.get("/visualization/{session_id}")
async def get_visualization_data(session_id: str, color_by: str = "scfocus_pseudotime"):
    """获取可视化数据"""
    return scfocus_service.get_visualization_data(session_id, color_by)


@router.get("/branches/{session_id}")
async def get_branch_data(session_id: str):
    """获取分支数据"""
    return scfocus_service.get_branch_data(session_id)
