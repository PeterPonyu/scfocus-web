"""结果获取路由"""
from fastapi import APIRouter
from fastapi.responses import FileResponse
from services.scfocus_service import scfocus_service
import os

router = APIRouter()

DATA_DIR = os.environ.get("DATA_DIR", "/tmp/scfocus_data")


@router.get("/data/{session_id}")
async def get_results(session_id: str):
    """获取分析结果数据"""
    return scfocus_service.get_results(session_id)


@router.get("/pseudotime/{session_id}")
async def get_pseudotime(session_id: str):
    """获取pseudotime数据"""
    result = scfocus_service.get_results(session_id)
    if not result.get("success"):
        return result
    
    return {
        "success": True,
        "pseudotime": result["data"]["pseudotime"],
        "n_cells": len(result["data"]["pseudotime"])
    }


@router.get("/entropy/{session_id}")
async def get_entropy(session_id: str):
    """获取entropy数据"""
    result = scfocus_service.get_results(session_id)
    if not result.get("success"):
        return result
    
    return {
        "success": True,
        "entropy": result["data"]["entropy"],
        "n_cells": len(result["data"]["entropy"])
    }


@router.get("/training-metrics/{session_id}")
async def get_training_metrics(session_id: str):
    """获取训练过程指标"""
    result = scfocus_service.get_results(session_id)
    if not result.get("success"):
        return result
    
    return {
        "success": True,
        "rewards": result["data"].get("training_rewards", []),
        "errors": result["data"].get("training_errors", [])
    }


@router.post("/export/{session_id}")
async def export_results(session_id: str, filename: str = "scfocus_results.h5ad"):
    """导出分析结果为h5ad文件"""
    output_path = os.path.join(DATA_DIR, f"{session_id}_{filename}")
    result = scfocus_service.export_results(session_id, output_path)
    
    if result.get("success"):
        return FileResponse(
            output_path,
            filename=filename,
            media_type="application/octet-stream"
        )
    
    return result


@router.get("/summary/{session_id}")
async def get_summary(session_id: str):
    """获取分析结果摘要"""
    status = scfocus_service.get_status(session_id)
    results = scfocus_service.get_results(session_id)
    branches = scfocus_service.get_branch_data(session_id)
    
    return {
        "success": True,
        "status": status,
        "has_results": results.get("success", False),
        "n_branches": branches.get("n_branches", 0) if branches.get("success") else 0,
        "simulated": results.get("data", {}).get("simulated", False) if results.get("success") else False
    }
