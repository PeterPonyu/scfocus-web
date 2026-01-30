"""
scFocus Web API - FastAPI后端
提供单细胞数据分析的强化学习聚焦分析接口
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from routers import data, analysis, training, results

app = FastAPI(
    title="scFocus Web API",
    description="单细胞强化学习聚焦分析Web服务 - 使用SAC算法进行谱系分支识别",
    version="1.0.0"
)

# 配置CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 注册路由
app.include_router(data.router, prefix="/api/data", tags=["数据管理"])
app.include_router(analysis.router, prefix="/api/analysis", tags=["数据分析"])
app.include_router(training.router, prefix="/api/training", tags=["模型训练"])
app.include_router(results.router, prefix="/api/results", tags=["结果获取"])


@app.get("/")
async def root():
    return {
        "message": "scFocus Web API",
        "description": "单细胞强化学习聚焦分析服务",
        "version": "1.0.0",
        "algorithm": "Soft Actor-Critic (SAC)",
        "docs": "/docs"
    }


@app.get("/health")
async def health_check():
    return {"status": "healthy"}
