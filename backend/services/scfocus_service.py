"""
scFocus服务封装
提供scFocus核心功能的封装接口
"""
import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple
import scanpy as sc
import anndata as ad
import os
import uuid
import json
from datetime import datetime

# 数据存储目录
DATA_DIR = os.environ.get("DATA_DIR", "/tmp/scfocus_data")
os.makedirs(DATA_DIR, exist_ok=True)


class ScFocusService:
    """scFocus服务类，封装所有分析功能"""
    
    def __init__(self):
        self.datasets: Dict[str, ad.AnnData] = {}
        self.focus_models: Dict[str, Any] = {}
        self.training_status: Dict[str, Dict] = {}
        self.results: Dict[str, Dict] = {}
    
    def create_session(self) -> str:
        """创建新的分析会话"""
        session_id = str(uuid.uuid4())
        self.training_status[session_id] = {
            "status": "initialized",
            "progress": 0,
            "message": "会话已创建",
            "created_at": datetime.now().isoformat()
        }
        return session_id
    
    def load_data(self, session_id: str, file_path: str, file_type: str = "h5ad") -> Dict[str, Any]:
        """
        加载单细胞数据
        
        Parameters
        ----------
        session_id : str
            会话ID
        file_path : str
            文件路径
        file_type : str
            文件类型 (h5ad, csv, mtx)
            
        Returns
        -------
        Dict
            数据概览信息
        """
        try:
            if file_type == "h5ad":
                adata = sc.read_h5ad(file_path)
            elif file_type == "csv":
                df = pd.read_csv(file_path, index_col=0)
                adata = ad.AnnData(df)
            elif file_type == "mtx":
                adata = sc.read_10x_mtx(os.path.dirname(file_path))
            else:
                raise ValueError(f"不支持的文件类型: {file_type}")
            
            self.datasets[session_id] = adata
            
            return {
                "success": True,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "obs_columns": list(adata.obs.columns),
                "var_columns": list(adata.var.columns),
                "obsm_keys": list(adata.obsm.keys()) if adata.obsm else [],
                "message": f"成功加载数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因"
            }
        except Exception as e:
            return {
                "success": False,
                "message": f"数据加载失败: {str(e)}"
            }
    
    def load_demo_data(self, session_id: str, dataset_id: str) -> Dict[str, Any]:
        """
        加载示例数据集
        
        Parameters
        ----------
        session_id : str
            会话ID
        dataset_id : str
            数据集ID: paul15, simulation_bifurcation, simulation_tree, pbmc3k
            
        Returns
        -------
        Dict
            数据概览信息
        """
        try:
            if dataset_id == "paul15":
                # 加载 Paul et al. 2015 小鼠造血数据
                adata = sc.datasets.paul15()
                self._update_status(session_id, "loading", 50, "加载 Paul15 小鼠造血数据...")
                
            elif dataset_id == "pbmc3k":
                # 加载 PBMC 3k 数据
                adata = sc.datasets.pbmc3k()
                self._update_status(session_id, "loading", 50, "加载 PBMC 3k 数据...")
                
            elif dataset_id == "simulation_bifurcation":
                # 生成分叉轨迹模拟数据
                adata = self._generate_bifurcation_simulation()
                self._update_status(session_id, "loading", 50, "生成分叉轨迹模拟数据...")
                
            elif dataset_id == "simulation_tree":
                # 生成树状轨迹模拟数据
                adata = self._generate_tree_simulation()
                self._update_status(session_id, "loading", 50, "生成树状轨迹模拟数据...")
                
            else:
                return {"success": False, "message": f"未知的数据集: {dataset_id}"}
            
            self.datasets[session_id] = adata
            self._update_status(session_id, "loaded", 100, f"数据集 {dataset_id} 加载完成")
            
            return {
                "success": True,
                "dataset_id": dataset_id,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "obs_columns": list(adata.obs.columns),
                "var_columns": list(adata.var.columns),
                "obsm_keys": list(adata.obsm.keys()) if adata.obsm else [],
                "message": f"成功加载示例数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因"
            }
        except Exception as e:
            self._update_status(session_id, "error", 0, f"加载失败: {str(e)}")
            return {"success": False, "message": f"数据加载失败: {str(e)}"}
    
    def _generate_bifurcation_simulation(self, n_cells: int = 1000, n_genes: int = 500) -> ad.AnnData:
        """生成分叉轨迹模拟数据"""
        np.random.seed(42)
        
        # 生成主干 + 两个分支
        n_trunk = n_cells // 3
        n_branch = (n_cells - n_trunk) // 2
        
        # 主干细胞 (早期)
        trunk_cells = np.random.randn(n_trunk, n_genes) * 0.5
        trunk_pseudotime = np.linspace(0, 0.4, n_trunk)
        trunk_branch = np.zeros(n_trunk)
        
        # 分支1细胞
        branch1_cells = np.random.randn(n_branch, n_genes) * 0.5 + np.array([1] * (n_genes // 2) + [-1] * (n_genes - n_genes // 2))
        branch1_pseudotime = np.linspace(0.4, 1.0, n_branch)
        branch1_branch = np.ones(n_branch)
        
        # 分支2细胞
        branch2_cells = np.random.randn(n_branch, n_genes) * 0.5 + np.array([-1] * (n_genes // 2) + [1] * (n_genes - n_genes // 2))
        branch2_pseudotime = np.linspace(0.4, 1.0, n_branch)
        branch2_branch = np.ones(n_branch) * 2
        
        # 合并
        X = np.vstack([trunk_cells, branch1_cells, branch2_cells])
        pseudotime = np.concatenate([trunk_pseudotime, branch1_pseudotime, branch2_pseudotime])
        branch = np.concatenate([trunk_branch, branch1_branch, branch2_branch])
        
        adata = ad.AnnData(X)
        adata.obs['pseudotime'] = pseudotime
        adata.obs['branch'] = branch.astype(int).astype(str)
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
        
        return adata
    
    def _generate_tree_simulation(self, n_cells: int = 2000, n_genes: int = 500) -> ad.AnnData:
        """生成树状轨迹模拟数据"""
        np.random.seed(42)
        
        # 生成具有多个分支的树状结构
        n_root = n_cells // 5
        n_branch1 = n_cells // 5
        n_branch2 = n_cells // 5
        n_leaf1 = n_cells // 5
        n_leaf2 = n_cells - n_root - n_branch1 - n_branch2 - n_leaf1
        
        # 根节点
        root_X = np.random.randn(n_root, n_genes) * 0.3
        root_pt = np.linspace(0, 0.25, n_root)
        root_branch = np.zeros(n_root)
        
        # 分支1
        branch1_X = np.random.randn(n_branch1, n_genes) * 0.3 + np.random.randn(n_genes) * 0.5
        branch1_pt = np.linspace(0.25, 0.5, n_branch1)
        branch1_branch = np.ones(n_branch1)
        
        # 分支2
        branch2_X = np.random.randn(n_branch2, n_genes) * 0.3 - np.random.randn(n_genes) * 0.5
        branch2_pt = np.linspace(0.25, 0.5, n_branch2)
        branch2_branch = np.ones(n_branch2) * 2
        
        # 叶子1 (从分支1延伸)
        leaf1_X = np.random.randn(n_leaf1, n_genes) * 0.3 + np.random.randn(n_genes) * 0.8
        leaf1_pt = np.linspace(0.5, 1.0, n_leaf1)
        leaf1_branch = np.ones(n_leaf1) * 3
        
        # 叶子2 (从分支2延伸)
        leaf2_X = np.random.randn(n_leaf2, n_genes) * 0.3 - np.random.randn(n_genes) * 0.8
        leaf2_pt = np.linspace(0.5, 1.0, n_leaf2)
        leaf2_branch = np.ones(n_leaf2) * 4
        
        # 合并
        X = np.vstack([root_X, branch1_X, branch2_X, leaf1_X, leaf2_X])
        pseudotime = np.concatenate([root_pt, branch1_pt, branch2_pt, leaf1_pt, leaf2_pt])
        branch = np.concatenate([root_branch, branch1_branch, branch2_branch, leaf1_branch, leaf2_branch])
        
        adata = ad.AnnData(X)
        adata.obs['pseudotime'] = pseudotime
        adata.obs['branch'] = branch.astype(int).astype(str)
        adata.obs['cell_type'] = pd.Categorical(['Root'] * n_root + ['Branch1'] * n_branch1 + 
                                                 ['Branch2'] * n_branch2 + ['Leaf1'] * n_leaf1 + ['Leaf2'] * n_leaf2)
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
        
        return adata
    
    def preprocess_data(
        self, 
        session_id: str,
        min_genes: int = 200,
        min_cells: int = 3,
        n_top_genes: int = 2000,
        n_pcs: int = 50,
        n_neighbors: int = 15,
        use_existing_embedding: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        预处理单细胞数据
        
        Parameters
        ----------
        session_id : str
            会话ID
        min_genes : int
            最小基因数阈值
        min_cells : int
            最小细胞数阈值
        n_top_genes : int
            高变基因数量
        n_pcs : int
            PCA维数
        n_neighbors : int
            邻居数量
        use_existing_embedding : str, optional
            使用现有embedding的key
            
        Returns
        -------
        Dict
            预处理结果信息
        """
        try:
            if session_id not in self.datasets:
                return {"success": False, "message": "未找到数据，请先加载数据"}
            
            adata = self.datasets[session_id]
            
            self._update_status(session_id, "preprocessing", 10, "开始预处理...")
            
            if use_existing_embedding and use_existing_embedding in adata.obsm:
                # 使用现有的embedding
                self._update_status(session_id, "preprocessing", 100, f"使用现有embedding: {use_existing_embedding}")
                return {
                    "success": True,
                    "n_cells": adata.n_obs,
                    "embedding_dim": adata.obsm[use_existing_embedding].shape[1],
                    "embedding_key": use_existing_embedding,
                    "message": f"使用现有embedding: {use_existing_embedding}"
                }
            
            # 质量控制
            self._update_status(session_id, "preprocessing", 20, "质量控制...")
            sc.pp.filter_cells(adata, min_genes=min_genes)
            sc.pp.filter_genes(adata, min_cells=min_cells)
            
            # 归一化
            self._update_status(session_id, "preprocessing", 40, "归一化...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            # 高变基因
            self._update_status(session_id, "preprocessing", 60, "识别高变基因...")
            sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars))
            adata_hvg = adata[:, adata.var.highly_variable].copy()
            
            # PCA
            self._update_status(session_id, "preprocessing", 80, "PCA降维...")
            sc.pp.scale(adata_hvg)
            sc.tl.pca(adata_hvg, n_comps=min(n_pcs, adata_hvg.n_vars - 1))
            
            # 邻居图和UMAP
            self._update_status(session_id, "preprocessing", 90, "计算UMAP...")
            sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors)
            sc.tl.umap(adata_hvg)
            
            # 更新数据
            adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
            adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
            self.datasets[session_id] = adata
            
            self._update_status(session_id, "preprocessed", 100, "预处理完成")
            
            return {
                "success": True,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "n_hvg": int(adata.var.highly_variable.sum()),
                "n_pcs": adata.obsm['X_pca'].shape[1],
                "has_umap": "X_umap" in adata.obsm,
                "message": "预处理完成"
            }
        except Exception as e:
            self._update_status(session_id, "error", 0, f"预处理失败: {str(e)}")
            return {"success": False, "message": f"预处理失败: {str(e)}"}
    
    def run_focus(
        self,
        session_id: str,
        embedding_key: str = "X_pca",
        hidden_dim: int = 128,
        n_agents: int = 8,
        max_steps: int = 5,
        pct_samples: float = 0.125,
        n_states: int = 2,
        num_episodes: int = 1000,
        batch_size: int = 64,
        meta_iterations: int = 3,
        resolution: float = 0.05
    ) -> Dict[str, Any]:
        """
        运行scFocus分析
        
        Parameters
        ----------
        session_id : str
            会话ID
        embedding_key : str
            用于分析的embedding key
        hidden_dim : int
            神经网络隐藏层维度
        n_agents : int
            并行环境/agent数量
        max_steps : int
            每episode最大步数
        pct_samples : float
            采样比例
        n_states : int
            状态变量数
        num_episodes : int
            训练episode数
        batch_size : int
            批次大小
        meta_iterations : int
            meta focusing迭代次数
        resolution : float
            合并分辨率参数
            
        Returns
        -------
        Dict
            分析结果
        """
        try:
            if session_id not in self.datasets:
                return {"success": False, "message": "未找到数据，请先加载数据"}
            
            adata = self.datasets[session_id]
            
            if embedding_key not in adata.obsm:
                return {"success": False, "message": f"未找到embedding: {embedding_key}"}
            
            # 获取latent space
            f = adata.obsm[embedding_key]
            
            self._update_status(session_id, "training", 5, "初始化scFocus模型...")
            
            try:
                # 尝试导入scFocus
                import sys
                sys.path.insert(0, os.path.expanduser("~/Projects/scFocus-source"))
                from scfocus import focus
                
                # 创建focus对象
                focus_model = focus(
                    f=f,
                    hidden_dim=hidden_dim,
                    n=n_agents,
                    max_steps=max_steps,
                    pct_samples=pct_samples,
                    n_states=n_states,
                    num_episodes=int(num_episodes),
                    batch_size=batch_size,
                    res=resolution
                )
                
                self._update_status(session_id, "training", 20, f"开始meta focusing ({meta_iterations}次迭代)...")
                
                # 运行meta focusing
                focus_model.meta_focusing(meta_iterations)
                
                self._update_status(session_id, "training", 80, "合并focus patterns...")
                
                # 合并focus patterns
                focus_model.merge_fp2()
                
                self._update_status(session_id, "training", 90, "计算pseudotime...")
                
                # 计算pseudotime
                focus_model.focus_diff()
                
                # 保存结果
                self.focus_models[session_id] = focus_model
                
                # 将结果保存到adata
                adata.obs['scfocus_pseudotime'] = focus_model.pseudotime
                adata.obs['scfocus_entropy'] = focus_model.entropy
                
                # 保存focus patterns
                for i, mfp in enumerate(focus_model.mfp):
                    for j in range(mfp.shape[1]):
                        adata.obs[f'scfocus_branch_{i}_{j}'] = mfp[:, j]
                
                self.datasets[session_id] = adata
                
                self._update_status(session_id, "completed", 100, "scFocus分析完成")
                
                # 计算结果统计
                n_branches = sum(mfp.shape[1] for mfp in focus_model.mfp)
                
                self.results[session_id] = {
                    "pseudotime": focus_model.pseudotime.tolist(),
                    "entropy": focus_model.entropy.tolist(),
                    "n_branches": n_branches,
                    "training_rewards": [r.tolist() for r in focus_model.r],
                    "training_errors": [e.tolist() for e in focus_model.e]
                }
                
                return {
                    "success": True,
                    "n_branches": n_branches,
                    "pseudotime_range": [float(focus_model.pseudotime.min()), float(focus_model.pseudotime.max())],
                    "entropy_range": [float(focus_model.entropy.min()), float(focus_model.entropy.max())],
                    "message": f"scFocus分析完成，识别到 {n_branches} 个谱系分支"
                }
                
            except ImportError as e:
                # 如果无法导入scFocus，使用模拟数据
                self._update_status(session_id, "training", 50, "使用模拟模式（未安装scFocus）...")
                
                n_cells = adata.n_obs
                
                # 生成模拟结果
                np.random.seed(42)
                pseudotime = np.random.rand(n_cells)
                entropy = np.random.rand(n_cells) * 2
                n_branches = np.random.randint(2, 5)
                
                adata.obs['scfocus_pseudotime'] = pseudotime
                adata.obs['scfocus_entropy'] = entropy
                
                for i in range(n_branches):
                    adata.obs[f'scfocus_branch_0_{i}'] = np.random.rand(n_cells)
                
                self.datasets[session_id] = adata
                
                self.results[session_id] = {
                    "pseudotime": pseudotime.tolist(),
                    "entropy": entropy.tolist(),
                    "n_branches": n_branches,
                    "training_rewards": [[float(np.random.rand())] * 100],
                    "training_errors": [[float(np.random.rand())] * 100],
                    "simulated": True
                }
                
                self._update_status(session_id, "completed", 100, "scFocus分析完成（模拟模式）")
                
                return {
                    "success": True,
                    "n_branches": n_branches,
                    "pseudotime_range": [0.0, 1.0],
                    "entropy_range": [0.0, 2.0],
                    "simulated": True,
                    "message": f"scFocus分析完成（模拟模式），识别到 {n_branches} 个谱系分支"
                }
                
        except Exception as e:
            self._update_status(session_id, "error", 0, f"分析失败: {str(e)}")
            return {"success": False, "message": f"分析失败: {str(e)}"}
    
    def get_results(self, session_id: str) -> Dict[str, Any]:
        """获取分析结果"""
        if session_id not in self.results:
            return {"success": False, "message": "未找到分析结果"}
        
        return {
            "success": True,
            "data": self.results[session_id]
        }
    
    def get_visualization_data(self, session_id: str, color_by: str = "scfocus_pseudotime") -> Dict[str, Any]:
        """
        获取可视化数据
        
        Parameters
        ----------
        session_id : str
            会话ID
        color_by : str
            着色依据的列名
            
        Returns
        -------
        Dict
            可视化数据
        """
        try:
            if session_id not in self.datasets:
                return {"success": False, "message": "未找到数据"}
            
            adata = self.datasets[session_id]
            
            if "X_umap" not in adata.obsm:
                return {"success": False, "message": "未找到UMAP坐标，请先进行预处理"}
            
            umap = adata.obsm["X_umap"]
            
            result = {
                "success": True,
                "x": umap[:, 0].tolist(),
                "y": umap[:, 1].tolist(),
                "cell_ids": adata.obs_names.tolist()
            }
            
            # 添加着色数据
            if color_by in adata.obs.columns:
                values = adata.obs[color_by].values
                if hasattr(values, 'tolist'):
                    result["colors"] = values.tolist()
                else:
                    result["colors"] = list(values)
                result["color_by"] = color_by
            
            # 添加可用的着色选项
            result["available_color_options"] = [
                col for col in adata.obs.columns 
                if col.startswith("scfocus_") or col in ["leiden", "louvain", "cell_type"]
            ]
            
            return result
        except Exception as e:
            return {"success": False, "message": f"获取可视化数据失败: {str(e)}"}
    
    def get_status(self, session_id: str) -> Dict[str, Any]:
        """获取会话状态"""
        if session_id not in self.training_status:
            return {"success": False, "message": "未找到会话"}
        
        return {
            "success": True,
            **self.training_status[session_id]
        }
    
    def get_branch_data(self, session_id: str) -> Dict[str, Any]:
        """获取分支数据"""
        try:
            if session_id not in self.datasets:
                return {"success": False, "message": "未找到数据"}
            
            adata = self.datasets[session_id]
            
            branch_columns = [col for col in adata.obs.columns if col.startswith("scfocus_branch_")]
            
            if not branch_columns:
                return {"success": False, "message": "未找到分支数据"}
            
            branches = {}
            for col in branch_columns:
                branches[col] = adata.obs[col].values.tolist()
            
            return {
                "success": True,
                "branches": branches,
                "n_branches": len(branch_columns)
            }
        except Exception as e:
            return {"success": False, "message": f"获取分支数据失败: {str(e)}"}
    
    def export_results(self, session_id: str, output_path: str) -> Dict[str, Any]:
        """导出分析结果"""
        try:
            if session_id not in self.datasets:
                return {"success": False, "message": "未找到数据"}
            
            adata = self.datasets[session_id]
            
            # 导出h5ad
            adata.write_h5ad(output_path)
            
            return {
                "success": True,
                "message": f"结果已导出到: {output_path}"
            }
        except Exception as e:
            return {"success": False, "message": f"导出失败: {str(e)}"}
    
    def _update_status(self, session_id: str, status: str, progress: int, message: str):
        """更新会话状态"""
        self.training_status[session_id] = {
            "status": status,
            "progress": progress,
            "message": message,
            "updated_at": datetime.now().isoformat()
        }


# 全局服务实例
scfocus_service = ScFocusService()
