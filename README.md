# scFocus Web

A web application for single-cell data analysis using the Soft Actor-Critic (SAC) reinforcement learning algorithm.

**🌐 Live Demo**: [https://scfocus-web.vercel.app](https://scfocus-web.vercel.app)

## Features

- **Single-cell Data Analysis**: Support for h5ad and csv format single-cell data upload and processing
- **SAC Reinforcement Learning**: Lineage branch identification using Soft Actor-Critic algorithm
- **Meta Focusing**: Multi-round iterative focusing analysis without prior knowledge
- **Interactive Visualization**: UMAP dimensionality reduction with multiple coloring options
- **Pseudotime Analysis**: Automatic cell developmental pseudotime calculation
- **Result Export**: Export analysis results in h5ad format
- **Bilingual Support**: English and Chinese interface

## Quick Start

### Online Usage

Visit [https://scfocus-web.vercel.app](https://scfocus-web.vercel.app) to use the application directly.

### Local Development

#### Backend

```bash
cd backend
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

#### Frontend

```bash
cd frontend
npm install
NEXT_PUBLIC_API_URL=http://localhost:8000 npm run dev
```

Visit http://localhost:3000 to access the application.

## Tech Stack

### Backend
- **FastAPI** - High-performance Python web framework
- **scFocus** - Single-cell reinforcement learning focusing analysis library
- **Scanpy/AnnData** - Single-cell data processing
- **PyTorch** - Deep learning framework

### Frontend
- **Next.js 14** - React full-stack framework
- **TailwindCSS** - CSS framework
- **React Query** - Data fetching
- **Recharts** - Chart library

## Usage Workflow

1. **Upload Data**: Upload single-cell data file (h5ad/csv) or load demo dataset
2. **Preprocess**: Set QC parameters, perform data preprocessing and dimensionality reduction
3. **Focus Analysis**: Configure SAC parameters, run meta focusing analysis
4. **View Results**: Visualize UMAP plots, view pseudotime and branch information
5. **Export Results**: Download h5ad file containing analysis results

## Deployment

### Backend (Render)

- **URL**: https://scfocus-api.onrender.com
- **Build Command**: `pip install -r requirements.txt`
- **Start Command**: `uvicorn main:app --host 0.0.0.0 --port $PORT`

### Frontend (Vercel)

- **URL**: https://scfocus-web.vercel.app
- **Environment Variable**: `NEXT_PUBLIC_API_URL=https://scfocus-api.onrender.com`

## API Documentation

Visit https://scfocus-api.onrender.com/docs for the complete API documentation.

### Main Endpoints

| Endpoint | Description |
|----------|-------------|
| `POST /api/data/session` | Create analysis session |
| `POST /api/data/upload/{session_id}` | Upload data file |
| `POST /api/analysis/preprocess/{session_id}` | Data preprocessing |
| `POST /api/training/focus/{session_id}` | Run scFocus analysis |
| `GET /api/results/data/{session_id}` | Get analysis results |

## Parameters

### Preprocessing Parameters
| Parameter | Description |
|-----------|-------------|
| `min_genes` | Minimum genes threshold |
| `min_cells` | Minimum cells threshold |
| `n_top_genes` | Number of highly variable genes |
| `n_pcs` | Number of PCA components |
| `n_neighbors` | Number of UMAP neighbors |

### scFocus Parameters
| Parameter | Description |
|-----------|-------------|
| `hidden_dim` | Neural network hidden layer dimension |
| `n_agents` | Number of parallel agents |
| `max_steps` | Maximum steps per episode |
| `pct_samples` | Sampling ratio |
| `num_episodes` | Number of training episodes |
| `meta_iterations` | Meta focusing iterations |

## License

MIT License

## Contact

- **Author**: Zeyu Fu
- **Email**: fuzeyu99@126.com
- `resolution`: 分支合并分辨率

## 许可证

MIT License

## 致谢

- scFocus原始算法: [PeterPonyu/scFocus](https://github.com/PeterPonyu/scFocus)
