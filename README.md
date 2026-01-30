# scFocus Web

基于Soft Actor-Critic强化学习算法的单细胞数据聚焦分析Web应用。

## 功能特性

- 🧬 **单细胞数据分析**: 支持h5ad、csv等格式的单细胞数据上传和处理
- 🤖 **SAC强化学习**: 使用Soft Actor-Critic算法进行谱系分支识别
- 🔬 **Meta Focusing**: 多轮迭代聚焦分析，无需先验知识
- 📊 **交互式可视化**: UMAP降维可视化，支持多种着色方式
- ⏱️ **Pseudotime分析**: 自动计算细胞发育pseudotime
- 📥 **结果导出**: 支持导出分析结果为h5ad格式

## 技术栈

### 后端
- FastAPI - 高性能Python Web框架
- scFocus - 单细胞强化学习聚焦分析库
- Scanpy/AnnData - 单细胞数据处理
- PyTorch - 深度学习框架

### 前端
- Next.js 14 - React全栈框架
- TailwindCSS - CSS框架
- React Query - 数据获取
- Recharts - 图表库

## 本地开发

### 后端

```bash
cd backend
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

### 前端

```bash
cd frontend
npm install
npm run dev
```

访问 http://localhost:3000 查看应用。

## 部署

### 后端部署 (Render)

1. 在Render创建新的Web Service
2. 连接GitHub仓库
3. 设置构建命令: `pip install -r requirements.txt`
4. 设置启动命令: `uvicorn main:app --host 0.0.0.0 --port $PORT`

### 前端部署 (Vercel)

1. 在Vercel导入项目
2. 设置环境变量 `NEXT_PUBLIC_API_URL` 为后端地址
3. 自动部署

## API文档

启动后端后访问 http://localhost:8000/docs 查看完整API文档。

### 主要接口

- `POST /api/data/session` - 创建分析会话
- `POST /api/data/upload/{session_id}` - 上传数据文件
- `POST /api/analysis/preprocess/{session_id}` - 数据预处理
- `POST /api/training/focus/{session_id}` - 运行scFocus分析
- `GET /api/results/data/{session_id}` - 获取分析结果

## 使用流程

1. **数据上传**: 上传单细胞数据文件(h5ad/csv)
2. **预处理**: 设置质控参数，进行数据预处理和降维
3. **Focus分析**: 配置SAC参数，运行meta focusing分析
4. **结果查看**: 可视化UMAP图，查看pseudotime和分支信息
5. **导出结果**: 下载包含分析结果的h5ad文件

## 参数说明

### 预处理参数
- `min_genes`: 最小基因数阈值
- `min_cells`: 最小细胞数阈值  
- `n_top_genes`: 高变基因数量
- `n_pcs`: PCA组分数
- `n_neighbors`: UMAP邻居数

### scFocus参数
- `hidden_dim`: 神经网络隐藏层维度
- `n_agents`: 并行Agent数量
- `max_steps`: 每episode最大步数
- `pct_samples`: 采样比例
- `num_episodes`: 训练轮数
- `meta_iterations`: Meta focusing迭代次数
- `resolution`: 分支合并分辨率

## 许可证

MIT License

## 致谢

- scFocus原始算法: [PeterPonyu/scFocus](https://github.com/PeterPonyu/scFocus)
