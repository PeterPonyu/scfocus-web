import axios from 'axios'

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'

const apiClient = axios.create({
  baseURL: API_URL,
  timeout: 300000, // 5分钟超时，因为训练可能很长
})

export const api = {
  // 会话管理
  createSession: async () => {
    const response = await apiClient.post('/api/data/session')
    return response.data
  },

  // 数据上传
  uploadFile: async (sessionId: string, file: File) => {
    const formData = new FormData()
    formData.append('file', file)
    formData.append('file_type', file.name.endsWith('.h5ad') ? 'h5ad' : 'csv')
    
    const response = await apiClient.post(`/api/data/upload/${sessionId}`, formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    })
    return response.data
  },

  // 加载演示数据
  loadDemoData: async (sessionId: string, dataset: string) => {
    const response = await apiClient.post(`/api/data/demo/${sessionId}`, { dataset })
    return response.data
  },

  // 获取演示数据列表
  getDemoDatasets: async () => {
    const response = await apiClient.get('/api/data/demo-datasets')
    return response.data
  },

  // 获取数据信息
  getDataInfo: async (sessionId: string) => {
    const response = await apiClient.get(`/api/data/info/${sessionId}`)
    return response.data
  },

  // 获取状态
  getStatus: async (sessionId: string) => {
    const response = await apiClient.get(`/api/data/status/${sessionId}`)
    return response.data
  },

  // 预处理
  preprocess: async (sessionId: string, params: any) => {
    const response = await apiClient.post(`/api/analysis/preprocess/${sessionId}`, params)
    return response.data
  },

  // 获取可视化数据
  getVisualization: async (sessionId: string, colorBy: string = 'scfocus_pseudotime') => {
    const response = await apiClient.get(`/api/analysis/visualization/${sessionId}`, {
      params: { color_by: colorBy }
    })
    return response.data
  },

  // 获取分支数据
  getBranches: async (sessionId: string) => {
    const response = await apiClient.get(`/api/analysis/branches/${sessionId}`)
    return response.data
  },

  // 运行Focus分析
  runFocus: async (sessionId: string, params: any) => {
    const response = await apiClient.post(`/api/training/focus/${sessionId}`, params)
    return response.data
  },

  // 获取训练状态
  getTrainingStatus: async (sessionId: string) => {
    const response = await apiClient.get(`/api/training/status/${sessionId}`)
    return response.data
  },

  // 获取结果
  getResults: async (sessionId: string) => {
    const response = await apiClient.get(`/api/results/data/${sessionId}`)
    return response.data
  },

  // 获取pseudotime
  getPseudotime: async (sessionId: string) => {
    const response = await apiClient.get(`/api/results/pseudotime/${sessionId}`)
    return response.data
  },

  // 获取entropy
  getEntropy: async (sessionId: string) => {
    const response = await apiClient.get(`/api/results/entropy/${sessionId}`)
    return response.data
  },

  // 获取训练指标
  getTrainingMetrics: async (sessionId: string) => {
    const response = await apiClient.get(`/api/results/training-metrics/${sessionId}`)
    return response.data
  },

  // 导出结果
  exportResults: async (sessionId: string, filename?: string) => {
    const response = await apiClient.post(`/api/results/export/${sessionId}`, null, {
      params: { filename },
      responseType: 'blob',
    })
    
    // 下载文件
    const url = window.URL.createObjectURL(new Blob([response.data]))
    const link = document.createElement('a')
    link.href = url
    link.setAttribute('download', filename || 'scfocus_results.h5ad')
    document.body.appendChild(link)
    link.click()
    link.remove()
    
    return { success: true }
  },

  // 获取摘要
  getSummary: async (sessionId: string) => {
    const response = await apiClient.get(`/api/results/summary/${sessionId}`)
    return response.data
  },
}
