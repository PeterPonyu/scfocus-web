'use client'

import { useState, useEffect } from 'react'
import { Upload, Settings, Play, BarChart3, Download, RefreshCw, Cpu, GitBranch, Activity, Zap } from 'lucide-react'
import FileUploader from '@/components/FileUploader'
import ScatterPlot from '@/components/ScatterPlot'
import PipelineSteps from '@/components/PipelineSteps'
import ParameterPanel from '@/components/ParameterPanel'
import ResultsPanel from '@/components/ResultsPanel'
import TrainingChart from '@/components/TrainingChart'
import { api } from '@/lib/api'

type Step = 'upload' | 'preprocess' | 'focus' | 'results'

interface Status {
  status: string
  progress: number
  message: string
}

export default function Home() {
  const [sessionId, setSessionId] = useState<string | null>(null)
  const [currentStep, setCurrentStep] = useState<Step>('upload')
  const [status, setStatus] = useState<Status | null>(null)
  const [dataInfo, setDataInfo] = useState<any>(null)
  const [visualizationData, setVisualizationData] = useState<any>(null)
  const [results, setResults] = useState<any>(null)
  const [colorBy, setColorBy] = useState<string>('scfocus_pseudotime')
  const [isLoading, setIsLoading] = useState(false)
  
  // 预处理参数
  const [preprocessParams, setPreprocessParams] = useState({
    min_genes: 200,
    min_cells: 3,
    n_top_genes: 2000,
    n_pcs: 50,
    n_neighbors: 15,
  })
  
  // Focus训练参数
  const [focusParams, setFocusParams] = useState({
    embedding_key: 'X_pca',
    hidden_dim: 128,
    n_agents: 8,
    max_steps: 5,
    pct_samples: 0.125,
    n_states: 2,
    num_episodes: 1000,
    batch_size: 64,
    meta_iterations: 3,
    resolution: 0.05,
  })

  // 初始化会话
  useEffect(() => {
    const initSession = async () => {
      try {
        const response = await api.createSession()
        setSessionId(response.session_id)
      } catch (error) {
        console.error('Failed to create session:', error)
      }
    }
    initSession()
  }, [])

  // 定时获取状态
  useEffect(() => {
    if (!sessionId || currentStep === 'upload' || currentStep === 'results') return
    
    const interval = setInterval(async () => {
      try {
        const statusResponse = await api.getStatus(sessionId)
        setStatus(statusResponse)
        
        if (statusResponse.status === 'completed') {
          setCurrentStep('results')
          await loadResults()
        }
      } catch (error) {
        console.error('Failed to get status:', error)
      }
    }, 2000)

    return () => clearInterval(interval)
  }, [sessionId, currentStep])

  const handleFileUpload = async (file: File) => {
    if (!sessionId) return
    
    setIsLoading(true)
    try {
      const response = await api.uploadFile(sessionId, file)
      if (response.success) {
        setDataInfo(response)
        setCurrentStep('preprocess')
      }
    } catch (error) {
      console.error('Upload failed:', error)
    } finally {
      setIsLoading(false)
    }
  }

  const handlePreprocess = async () => {
    if (!sessionId) return
    
    setIsLoading(true)
    setStatus({ status: 'preprocessing', progress: 0, message: '开始预处理...' })
    
    try {
      const response = await api.preprocess(sessionId, preprocessParams)
      if (response.success) {
        setDataInfo((prev: any) => ({ ...prev, ...response }))
        setCurrentStep('focus')
        
        // 获取可视化数据
        const vizData = await api.getVisualization(sessionId, 'scfocus_pseudotime')
        setVisualizationData(vizData)
      }
    } catch (error) {
      console.error('Preprocess failed:', error)
    } finally {
      setIsLoading(false)
      setStatus(null)
    }
  }

  const handleRunFocus = async () => {
    if (!sessionId) return
    
    setIsLoading(true)
    setStatus({ status: 'training', progress: 0, message: '初始化scFocus模型...' })
    
    try {
      const response = await api.runFocus(sessionId, focusParams)
      if (response.success) {
        setResults(response)
        await loadResults()
        setCurrentStep('results')
      }
    } catch (error) {
      console.error('Focus training failed:', error)
    } finally {
      setIsLoading(false)
    }
  }

  const loadResults = async () => {
    if (!sessionId) return
    
    try {
      const [vizData, resultsData, branchData] = await Promise.all([
        api.getVisualization(sessionId, colorBy),
        api.getResults(sessionId),
        api.getBranches(sessionId),
      ])
      
      setVisualizationData(vizData)
      setResults({ ...resultsData.data, branches: branchData })
    } catch (error) {
      console.error('Failed to load results:', error)
    }
  }

  const handleColorChange = async (newColorBy: string) => {
    if (!sessionId) return
    setColorBy(newColorBy)
    
    try {
      const vizData = await api.getVisualization(sessionId, newColorBy)
      setVisualizationData(vizData)
    } catch (error) {
      console.error('Failed to update visualization:', error)
    }
  }

  const handleExport = async () => {
    if (!sessionId) return
    
    try {
      await api.exportResults(sessionId)
    } catch (error) {
      console.error('Export failed:', error)
    }
  }

  const steps = [
    { id: 'upload', label: '数据上传', icon: Upload },
    { id: 'preprocess', label: '预处理', icon: Settings },
    { id: 'focus', label: 'Focus分析', icon: GitBranch },
    { id: 'results', label: '结果展示', icon: BarChart3 },
  ]

  return (
    <main className="min-h-screen p-6">
      {/* 头部 */}
      <header className="mb-8">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="p-3 bg-gradient-to-br from-green-500 to-purple-600 rounded-xl shadow-lg">
              <GitBranch className="w-8 h-8 text-white" />
            </div>
            <div>
              <h1 className="text-3xl font-bold text-gray-800">scFocus Web</h1>
              <p className="text-sm text-gray-500">单细胞强化学习聚焦分析平台 · SAC算法</p>
            </div>
          </div>
          
          <div className="flex items-center gap-4">
            <div className="flex items-center gap-2 px-4 py-2 bg-white rounded-lg shadow-sm">
              <Cpu className="w-4 h-4 text-green-500" />
              <span className="text-sm text-gray-600">
                会话: {sessionId ? sessionId.slice(0, 8) : '初始化中...'}
              </span>
            </div>
            
            {currentStep === 'results' && (
              <button
                onClick={handleExport}
                className="flex items-center gap-2 px-4 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-colors shadow-sm"
              >
                <Download className="w-4 h-4" />
                导出结果
              </button>
            )}
          </div>
        </div>
      </header>

      {/* 进度步骤 */}
      <PipelineSteps steps={steps} currentStep={currentStep} />

      {/* 状态栏 */}
      {status && status.status !== 'completed' && (
        <div className="mb-6 p-4 bg-white rounded-xl shadow-sm border border-gray-100 animate-fade-in">
          <div className="flex items-center gap-3 mb-2">
            <RefreshCw className="w-5 h-5 text-green-500 animate-spin" />
            <span className="font-medium text-gray-700">{status.message}</span>
          </div>
          <div className="w-full bg-gray-100 rounded-full h-2">
            <div
              className="bg-gradient-to-r from-green-500 to-purple-500 h-2 rounded-full transition-all duration-500"
              style={{ width: `${status.progress}%` }}
            />
          </div>
          <p className="text-sm text-gray-500 mt-1">{status.progress}% 完成</p>
        </div>
      )}

      {/* 主内容区 */}
      <div className="grid grid-cols-12 gap-6">
        {/* 左侧控制面板 */}
        <div className="col-span-4 space-y-6">
          {/* 数据上传 */}
          {currentStep === 'upload' && (
            <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
              <h2 className="text-lg font-semibold text-gray-800 mb-4 flex items-center gap-2">
                <Upload className="w-5 h-5 text-green-500" />
                上传单细胞数据
              </h2>
              <FileUploader onUpload={handleFileUpload} isLoading={isLoading} />
              <div className="mt-4 text-sm text-gray-500">
                <p>支持格式: .h5ad, .csv</p>
                <p>建议: 使用预处理后的AnnData文件效果更佳</p>
              </div>
            </div>
          )}

          {/* 预处理参数 */}
          {currentStep === 'preprocess' && (
            <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
              <h2 className="text-lg font-semibold text-gray-800 mb-4 flex items-center gap-2">
                <Settings className="w-5 h-5 text-green-500" />
                预处理参数
              </h2>
              <ParameterPanel
                params={preprocessParams}
                setParams={setPreprocessParams}
                paramConfig={[
                  { key: 'min_genes', label: '最小基因数', type: 'number', min: 0, max: 1000 },
                  { key: 'min_cells', label: '最小细胞数', type: 'number', min: 0, max: 100 },
                  { key: 'n_top_genes', label: '高变基因数', type: 'number', min: 500, max: 5000 },
                  { key: 'n_pcs', label: 'PCA维数', type: 'number', min: 10, max: 100 },
                  { key: 'n_neighbors', label: '邻居数', type: 'number', min: 5, max: 50 },
                ]}
              />
              <button
                onClick={handlePreprocess}
                disabled={isLoading}
                className="w-full mt-4 py-3 bg-gradient-to-r from-green-500 to-green-600 text-white rounded-lg font-medium hover:from-green-600 hover:to-green-700 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
              >
                {isLoading ? (
                  <RefreshCw className="w-5 h-5 animate-spin" />
                ) : (
                  <Play className="w-5 h-5" />
                )}
                {isLoading ? '处理中...' : '开始预处理'}
              </button>
            </div>
          )}

          {/* Focus训练参数 */}
          {currentStep === 'focus' && (
            <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
              <h2 className="text-lg font-semibold text-gray-800 mb-4 flex items-center gap-2">
                <GitBranch className="w-5 h-5 text-purple-500" />
                scFocus参数
              </h2>
              <ParameterPanel
                params={focusParams}
                setParams={setFocusParams}
                paramConfig={[
                  { key: 'hidden_dim', label: '隐藏层维度', type: 'number', min: 32, max: 512 },
                  { key: 'n_agents', label: 'Agent数量', type: 'number', min: 4, max: 16 },
                  { key: 'max_steps', label: '最大步数', type: 'number', min: 3, max: 10 },
                  { key: 'pct_samples', label: '采样比例', type: 'number', min: 0.05, max: 0.5, step: 0.025 },
                  { key: 'n_states', label: '状态变量数', type: 'number', min: 2, max: 10 },
                  { key: 'num_episodes', label: '训练轮数', type: 'number', min: 100, max: 5000 },
                  { key: 'batch_size', label: '批次大小', type: 'number', min: 16, max: 256 },
                  { key: 'meta_iterations', label: 'Meta迭代次数', type: 'number', min: 1, max: 10 },
                  { key: 'resolution', label: '合并分辨率', type: 'number', min: 0.01, max: 0.2, step: 0.01 },
                ]}
              />
              <button
                onClick={handleRunFocus}
                disabled={isLoading}
                className="w-full mt-4 py-3 bg-gradient-to-r from-purple-500 to-purple-600 text-white rounded-lg font-medium hover:from-purple-600 hover:to-purple-700 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
              >
                {isLoading ? (
                  <RefreshCw className="w-5 h-5 animate-spin" />
                ) : (
                  <Zap className="w-5 h-5" />
                )}
                {isLoading ? '训练中...' : '运行Focus分析'}
              </button>
            </div>
          )}

          {/* 数据信息 */}
          {dataInfo && (
            <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
              <h2 className="text-lg font-semibold text-gray-800 mb-4 flex items-center gap-2">
                <Activity className="w-5 h-5 text-blue-500" />
                数据信息
              </h2>
              <div className="space-y-3">
                <div className="flex justify-between items-center p-3 bg-gray-50 rounded-lg">
                  <span className="text-gray-600">细胞数</span>
                  <span className="font-semibold text-gray-800">{dataInfo.n_cells?.toLocaleString()}</span>
                </div>
                <div className="flex justify-between items-center p-3 bg-gray-50 rounded-lg">
                  <span className="text-gray-600">基因数</span>
                  <span className="font-semibold text-gray-800">{dataInfo.n_genes?.toLocaleString()}</span>
                </div>
                {dataInfo.n_hvg && (
                  <div className="flex justify-between items-center p-3 bg-gray-50 rounded-lg">
                    <span className="text-gray-600">高变基因</span>
                    <span className="font-semibold text-gray-800">{dataInfo.n_hvg?.toLocaleString()}</span>
                  </div>
                )}
                {dataInfo.obsm_keys && dataInfo.obsm_keys.length > 0 && (
                  <div className="p-3 bg-gray-50 rounded-lg">
                    <span className="text-gray-600">Embeddings:</span>
                    <div className="flex flex-wrap gap-1 mt-1">
                      {dataInfo.obsm_keys.map((key: string) => (
                        <span key={key} className="px-2 py-1 bg-green-100 text-green-700 text-xs rounded">
                          {key}
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </div>
          )}
        </div>

        {/* 右侧可视化区域 */}
        <div className="col-span-8 space-y-6">
          {/* 散点图可视化 */}
          {visualizationData && visualizationData.success && (
            <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
              <div className="flex items-center justify-between mb-4">
                <h2 className="text-lg font-semibold text-gray-800 flex items-center gap-2">
                  <BarChart3 className="w-5 h-5 text-green-500" />
                  UMAP可视化
                </h2>
                {visualizationData.available_color_options && (
                  <select
                    value={colorBy}
                    onChange={(e) => handleColorChange(e.target.value)}
                    className="px-3 py-2 bg-gray-50 border border-gray-200 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-green-500"
                  >
                    {visualizationData.available_color_options.map((option: string) => (
                      <option key={option} value={option}>
                        {option}
                      </option>
                    ))}
                  </select>
                )}
              </div>
              <ScatterPlot
                x={visualizationData.x}
                y={visualizationData.y}
                colors={visualizationData.colors}
                colorBy={colorBy}
              />
            </div>
          )}

          {/* 结果展示 */}
          {currentStep === 'results' && results && (
            <>
              <ResultsPanel results={results} />
              
              {results.training_rewards && (
                <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
                  <h2 className="text-lg font-semibold text-gray-800 mb-4 flex items-center gap-2">
                    <Activity className="w-5 h-5 text-purple-500" />
                    训练过程
                  </h2>
                  <TrainingChart
                    rewards={results.training_rewards}
                    errors={results.training_errors}
                  />
                </div>
              )}
            </>
          )}

          {/* 空状态提示 */}
          {currentStep === 'upload' && (
            <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-12 text-center animate-fade-in">
              <div className="inline-flex items-center justify-center w-20 h-20 bg-gradient-to-br from-green-100 to-purple-100 rounded-full mb-6">
                <GitBranch className="w-10 h-10 text-green-600" />
              </div>
              <h3 className="text-xl font-semibold text-gray-800 mb-2">开始scFocus分析</h3>
              <p className="text-gray-500 max-w-md mx-auto">
                上传您的单细胞数据，使用Soft Actor-Critic强化学习算法识别谱系分支并计算pseudotime
              </p>
              <div className="mt-6 flex justify-center gap-4">
                <div className="flex items-center gap-2 text-sm text-gray-500">
                  <div className="w-2 h-2 bg-green-500 rounded-full"></div>
                  SAC算法
                </div>
                <div className="flex items-center gap-2 text-sm text-gray-500">
                  <div className="w-2 h-2 bg-purple-500 rounded-full"></div>
                  Meta Focusing
                </div>
                <div className="flex items-center gap-2 text-sm text-gray-500">
                  <div className="w-2 h-2 bg-blue-500 rounded-full"></div>
                  无需先验知识
                </div>
              </div>
            </div>
          )}
        </div>
      </div>
    </main>
  )
}
