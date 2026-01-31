"""数据管理路由"""
from fastapi import APIRouter, UploadFile, File, Form, Query
from pydantic import BaseModel
from typing import Optional, List, Literal
import os
import shutil
from services.scfocus_service import scfocus_service

router = APIRouter()

DATA_DIR = os.environ.get("DATA_DIR", "/tmp/scfocus_data")
os.makedirs(DATA_DIR, exist_ok=True)


class CreateSessionResponse(BaseModel):
    session_id: str
    message: str


class DataLoadResponse(BaseModel):
    success: bool
    message: str
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    obs_columns: Optional[List[str]] = None
    var_columns: Optional[List[str]] = None
    obsm_keys: Optional[List[str]] = None


class DemoDatasetInfo(BaseModel):
    """示例数据集信息"""
    id: str
    name: str
    description: str
    n_cells: int
    n_genes: int
    has_pseudotime: bool
    source: str


@router.get("/demo-datasets", response_model=List[DemoDatasetInfo])
async def list_demo_datasets():
    """获取可用的示例数据集列表"""
    return [
        DemoDatasetInfo(
            id="paul15",
            name="Paul et al. 2015 (小鼠造血)",
            description="小鼠骨髓造血祖细胞分化轨迹数据，包含2730个细胞，常用于轨迹推断基准测试",
            n_cells=2730,
            n_genes=3451,
            has_pseudotime=True,
            source="scanpy.datasets.paul15()"
        ),
        DemoDatasetInfo(
            id="simulation_bifurcation",
            name="模拟数据 - 分叉轨迹",
            description="包含单个分叉点的模拟单细胞轨迹数据，适合测试谱系分支识别",
            n_cells=1000,
            n_genes=500,
            has_pseudotime=True,
            source="simulation"
        ),
        DemoDatasetInfo(
            id="simulation_tree",
            name="模拟数据 - 树状轨迹",
            description="包含多个分支的树状轨迹模拟数据，适合测试复杂谱系结构",
            n_cells=2000,
            n_genes=500,
            has_pseudotime=True,
            source="simulation"
        ),
        DemoDatasetInfo(
            id="pbmc3k",
            name="PBMC 3k (人类外周血)",
            description="10X Genomics提供的3000个人类外周血单核细胞数据",
            n_cells=2700,
            n_genes=13714,
            has_pseudotime=False,
            source="scanpy.datasets.pbmc3k()"
        ),
    ]


@router.post("/demo/{session_id}/{dataset_id}")
async def load_demo_data(
    session_id: str,
    dataset_id: Literal["paul15", "simulation_bifurcation", "simulation_tree", "pbmc3k"]
):
    """加载示例数据集"""
    result = scfocus_service.load_demo_data(session_id, dataset_id)
    return result


@router.post("/session", response_model=CreateSessionResponse)
async def create_session():
    """创建新的分析会话"""
    session_id = scfocus_service.create_session()
    return CreateSessionResponse(
        session_id=session_id,
        message="会话创建成功"
    )


@router.post("/upload/{session_id}")
async def upload_file(
    session_id: str,
    file: UploadFile = File(...),
    file_type: str = Form("h5ad")
):
    """上传数据文件"""
    # 保存文件
    file_path = os.path.join(DATA_DIR, f"{session_id}_{file.filename}")
    
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    
    # 加载数据
    result = scfocus_service.load_data(session_id, file_path, file_type)
    
    return result


@router.get("/info/{session_id}")
async def get_data_info(session_id: str):
    """获取数据信息"""
    if session_id not in scfocus_service.datasets:
        return {"success": False, "message": "未找到数据"}
    
    adata = scfocus_service.datasets[session_id]
    
    return {
        "success": True,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "obsm_keys": list(adata.obsm.keys()) if adata.obsm else []
    }


@router.get("/status/{session_id}")
async def get_status(session_id: str):
    """获取会话状态"""
    return scfocus_service.get_status(session_id)
