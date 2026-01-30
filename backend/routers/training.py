"""模型训练路由"""
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import Optional
from services.scfocus_service import scfocus_service

router = APIRouter()


class FocusTrainingRequest(BaseModel):
    embedding_key: str = "X_pca"
    hidden_dim: int = 128
    n_agents: int = 8
    max_steps: int = 5
    pct_samples: float = 0.125
    n_states: int = 2
    num_episodes: int = 1000
    batch_size: int = 64
    meta_iterations: int = 3
    resolution: float = 0.05


class TrainingResponse(BaseModel):
    success: bool
    message: str
    n_branches: Optional[int] = None
    simulated: Optional[bool] = None


@router.post("/focus/{session_id}", response_model=TrainingResponse)
async def run_focus_training(session_id: str, request: FocusTrainingRequest):
    """运行scFocus聚焦分析"""
    result = scfocus_service.run_focus(
        session_id=session_id,
        embedding_key=request.embedding_key,
        hidden_dim=request.hidden_dim,
        n_agents=request.n_agents,
        max_steps=request.max_steps,
        pct_samples=request.pct_samples,
        n_states=request.n_states,
        num_episodes=request.num_episodes,
        batch_size=request.batch_size,
        meta_iterations=request.meta_iterations,
        resolution=request.resolution
    )
    return result


@router.get("/status/{session_id}")
async def get_training_status(session_id: str):
    """获取训练状态"""
    return scfocus_service.get_status(session_id)


@router.post("/focus-async/{session_id}")
async def run_focus_training_async(
    session_id: str, 
    request: FocusTrainingRequest,
    background_tasks: BackgroundTasks
):
    """异步运行scFocus聚焦分析（用于长时间训练）"""
    
    def run_training():
        scfocus_service.run_focus(
            session_id=session_id,
            embedding_key=request.embedding_key,
            hidden_dim=request.hidden_dim,
            n_agents=request.n_agents,
            max_steps=request.max_steps,
            pct_samples=request.pct_samples,
            n_states=request.n_states,
            num_episodes=request.num_episodes,
            batch_size=request.batch_size,
            meta_iterations=request.meta_iterations,
            resolution=request.resolution
        )
    
    background_tasks.add_task(run_training)
    
    return {
        "success": True,
        "message": "训练任务已提交，请通过 /status 接口查询进度"
    }
