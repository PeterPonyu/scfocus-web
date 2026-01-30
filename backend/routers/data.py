"""数据管理路由"""
from fastapi import APIRouter, UploadFile, File, Form
from pydantic import BaseModel
from typing import Optional, List
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
