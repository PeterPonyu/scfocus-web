export const translations = {
  en: {
    // Header
    appName: 'scFocus Web',
    appDescription: 'Single-cell Reinforcement Learning Focusing Analysis Platform',
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
      dragDrop: 'Drag and drop files here, or',
      clickSelect: 'click to select',
      supportedFormats: 'Supported formats: .h5ad, .csv',
      suggestion: 'Preprocessed AnnData files are recommended for optimal results.',
      uploadButton: 'Upload Data',
      uploading: 'Uploading...',
    },
    
    // Preprocess
    preprocess: {
      title: 'Preprocessing Parameters',
      minGenes: 'Min Genes per Cell',
      minCells: 'Min Cells per Gene',
      topGenes: 'Number of HVGs',
      nPcs: 'PCA Components',
      nNeighbors: 'Number of Neighbors',
      startButton: 'Run Preprocessing',
      processing: 'Processing...',
    },
    
    // Focus
    focus: {
      title: 'scFocus Parameters',
      hiddenDim: 'Hidden Dimension',
      nAgents: 'Number of Agents',
      maxSteps: 'Max Steps',
      pctSamples: 'Sample Fraction',
      nStates: 'Number of States',
      numEpisodes: 'Training Episodes',
      batchSize: 'Batch Size',
      metaIterations: 'Meta Iterations',
      resolution: 'Merge Resolution',
      runButton: 'Run Analysis',
      training: 'Training...',
    },
    
    // Data Info
    dataInfo: {
      title: 'Data Summary',
      cells: 'Cells',
      genes: 'Genes',
      hvg: 'HVGs',
      embeddings: 'Embeddings',
    },
    
    // Visualization
    visualization: {
      title: 'UMAP Visualization',
    },
    
    // Results
    results: {
      title: 'Analysis Results',
      summary: 'Results Summary',
      branches: 'Identified Branches',
      pseudotimeRange: 'Pseudotime Range',
      entropyRange: 'Entropy Range',
      analysisMode: 'Analysis Mode',
      sacAlgorithm: 'SAC Algorithm',
      simulatedMode: 'Simulation Mode',
      simulatedWarning: 'Running in simulation mode (scFocus library not installed). Install scFocus on the backend for actual analysis.',
      trainingProcess: 'Training Progress',
    },
    
    // Welcome
    welcome: {
      title: 'scFocus Analysis',
      description: 'Upload single-cell data to identify lineage branches and compute pseudotime using the Soft Actor-Critic reinforcement learning algorithm.',
      features: {
        sac: 'SAC Algorithm',
        meta: 'Meta Focusing',
        noPrior: 'No Prior Required',
      },
    },
    
    // Status
    status: {
      starting: 'Starting preprocessing...',
      initializing: 'Initializing model...',
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
    appDescription: '单细胞强化学习聚焦分析平台',
    session: '会话',
    initializing: '初始化中...',
    exportResults: '导出结果',
    
    // Steps
    steps: {
      upload: '数据上传',
      preprocess: '预处理',
      focus: 'Focus分析',
      results: '结果',
    },
    
    // Upload
    upload: {
      title: '上传单细胞数据',
      dragDrop: '拖放文件到此处，或',
      clickSelect: '点击选择',
      supportedFormats: '支持格式: .h5ad, .csv',
      suggestion: '建议使用预处理后的AnnData文件以获得更好的分析效果。',
      uploadButton: '上传数据',
      uploading: '上传中...',
    },
    
    // Preprocess
    preprocess: {
      title: '预处理参数',
      minGenes: '每个细胞最小基因数',
      minCells: '每个基因最小细胞数',
      topGenes: '高变基因数量',
      nPcs: 'PCA主成分数',
      nNeighbors: '邻居数量',
      startButton: '运行预处理',
      processing: '处理中...',
    },
    
    // Focus
    focus: {
      title: 'scFocus参数',
      hiddenDim: '隐藏层维度',
      nAgents: 'Agent数量',
      maxSteps: '最大步数',
      pctSamples: '采样比例',
      nStates: '状态数量',
      numEpisodes: '训练轮数',
      batchSize: '批次大小',
      metaIterations: 'Meta迭代次数',
      resolution: '合并分辨率',
      runButton: '运行分析',
      training: '训练中...',
    },
    
    // Data Info
    dataInfo: {
      title: '数据概要',
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
      title: '分析结果',
      summary: '结果概要',
      branches: '识别的分支数',
      pseudotimeRange: 'Pseudotime范围',
      entropyRange: 'Entropy范围',
      analysisMode: '分析模式',
      sacAlgorithm: 'SAC算法',
      simulatedMode: '模拟模式',
      simulatedWarning: '当前为模拟模式（未安装scFocus库）。实际分析需要在后端安装scFocus。',
      trainingProcess: '训练进度',
    },
    
    // Welcome
    welcome: {
      title: 'scFocus分析',
      description: '上传单细胞数据，使用Soft Actor-Critic强化学习算法识别谱系分支并计算pseudotime。',
      features: {
        sac: 'SAC算法',
        meta: 'Meta Focusing',
        noPrior: '无需先验',
      },
    },
    
    // Status
    status: {
      starting: '开始预处理...',
      initializing: '初始化模型...',
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
