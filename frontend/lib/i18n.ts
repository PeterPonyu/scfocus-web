export const translations = {
  en: {
    // Header
    appName: 'scFocus Web',
    appDescription: 'Single-cell RL Focusing Analysis Platform · SAC Algorithm',
    session: 'Session',
    initializing: 'Initializing...',
    exportResults: 'Export Results',
    
    // Steps
    steps: {
      upload: 'Data Upload',
      preprocess: 'Preprocess',
      focus: 'Focus Analysis',
      results: 'Results',
    },
    
    // Upload
    upload: {
      title: 'Upload Single-cell Data',
      dragDrop: 'Drag and drop files here or',
      clickSelect: 'click to select',
      supportedFormats: 'Supported formats: .h5ad, .csv',
      suggestion: 'Suggestion: Use preprocessed AnnData files for better results',
      uploadButton: 'Upload and Load Data',
      uploading: 'Uploading...',
    },
    
    // Preprocess
    preprocess: {
      title: 'Preprocessing Parameters',
      minGenes: 'Min Genes',
      minCells: 'Min Cells',
      topGenes: 'HVG Count',
      nPcs: 'PCA Dimensions',
      nNeighbors: 'Neighbors',
      startButton: 'Start Preprocessing',
      processing: 'Processing...',
    },
    
    // Focus
    focus: {
      title: 'scFocus Parameters',
      hiddenDim: 'Hidden Dimensions',
      nAgents: 'Agent Count',
      maxSteps: 'Max Steps',
      pctSamples: 'Sample Ratio',
      nStates: 'State Variables',
      numEpisodes: 'Training Episodes',
      batchSize: 'Batch Size',
      metaIterations: 'Meta Iterations',
      resolution: 'Merge Resolution',
      runButton: 'Run Focus Analysis',
      training: 'Training...',
    },
    
    // Data Info
    dataInfo: {
      title: 'Data Info',
      cells: 'Cells',
      genes: 'Genes',
      hvg: 'HVG',
      embeddings: 'Embeddings',
    },
    
    // Visualization
    visualization: {
      title: 'UMAP Visualization',
    },
    
    // Results
    results: {
      title: 'Analysis Results Summary',
      summary: 'Analysis Results Summary',
      branches: 'Identified Branches',
      pseudotimeRange: 'Pseudotime Range',
      entropyRange: 'Entropy Range',
      analysisMode: 'Analysis Mode',
      sacAlgorithm: 'SAC Algorithm',
      simulatedMode: 'Simulated Mode',
      simulatedWarning: '⚠️ Currently using simulated data mode (scFocus library not installed). Actual analysis requires scFocus installation on the backend.',
      trainingProcess: 'Training Process',
    },
    
    // Welcome
    welcome: {
      title: 'Start scFocus Analysis',
      description: 'Upload your single-cell data, use Soft Actor-Critic reinforcement learning algorithm to identify lineage branches and calculate pseudotime',
      features: {
        sac: 'SAC Algorithm',
        meta: 'Meta Focusing',
        noPrior: 'No Prior Knowledge Required',
      },
    },
    
    // Status
    status: {
      starting: 'Starting preprocessing...',
      initializing: 'Initializing scFocus model...',
      complete: 'Complete',
    },
    
    // Language
    language: 'Language',
    chinese: '中文',
    english: 'English',
  },
  
  zh: {
    // Header
    appName: 'scFocus Web',
    appDescription: '单细胞强化学习聚焦分析平台 · SAC算法',
    session: '会话',
    initializing: '初始化中...',
    exportResults: '导出结果',
    
    // Steps
    steps: {
      upload: '数据上传',
      preprocess: '预处理',
      focus: 'Focus分析',
      results: '结果展示',
    },
    
    // Upload
    upload: {
      title: '上传单细胞数据',
      dragDrop: '拖放文件到此处或',
      clickSelect: '点击选择',
      supportedFormats: '支持格式: .h5ad, .csv',
      suggestion: '建议: 使用预处理后的AnnData文件效果更佳',
      uploadButton: '上传并加载数据',
      uploading: '上传中...',
    },
    
    // Preprocess
    preprocess: {
      title: '预处理参数',
      minGenes: '最小基因数',
      minCells: '最小细胞数',
      topGenes: '高变基因数',
      nPcs: 'PCA维数',
      nNeighbors: '邻居数',
      startButton: '开始预处理',
      processing: '处理中...',
    },
    
    // Focus
    focus: {
      title: 'scFocus参数',
      hiddenDim: '隐藏层维度',
      nAgents: 'Agent数量',
      maxSteps: '最大步数',
      pctSamples: '采样比例',
      nStates: '状态变量数',
      numEpisodes: '训练轮数',
      batchSize: '批次大小',
      metaIterations: 'Meta迭代次数',
      resolution: '合并分辨率',
      runButton: '运行Focus分析',
      training: '训练中...',
    },
    
    // Data Info
    dataInfo: {
      title: '数据信息',
      cells: '细胞数',
      genes: '基因数',
      hvg: '高变基因',
      embeddings: 'Embeddings',
    },
    
    // Visualization
    visualization: {
      title: 'UMAP可视化',
    },
    
    // Results
    results: {
      title: '分析结果摘要',
      summary: '分析结果摘要',
      branches: '识别分支数',
      pseudotimeRange: 'Pseudotime范围',
      entropyRange: 'Entropy范围',
      analysisMode: '分析模式',
      sacAlgorithm: 'SAC算法',
      simulatedMode: '模拟模式',
      simulatedWarning: '⚠️ 当前使用模拟数据模式（未安装scFocus库）。实际分析需要在后端环境中安装scFocus。',
      trainingProcess: '训练过程',
    },
    
    // Welcome
    welcome: {
      title: '开始scFocus分析',
      description: '上传您的单细胞数据，使用Soft Actor-Critic强化学习算法识别谱系分支并计算pseudotime',
      features: {
        sac: 'SAC算法',
        meta: 'Meta Focusing',
        noPrior: '无需先验知识',
      },
    },
    
    // Status
    status: {
      starting: '开始预处理...',
      initializing: '初始化scFocus模型...',
      complete: '完成',
    },
    
    // Language
    language: '语言',
    chinese: '中文',
    english: 'English',
  },
}

export type Language = 'en' | 'zh'
export type TranslationKey = keyof typeof translations.en
