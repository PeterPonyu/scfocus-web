'use client'

import { GitBranch, Clock, Activity, Sparkles } from 'lucide-react'

interface ResultsPanelProps {
  results: any
}

export default function ResultsPanel({ results }: ResultsPanelProps) {
  if (!results) return null

  const stats = [
    {
      label: '识别分支数',
      value: results.n_branches || results.branches?.n_branches || 0,
      icon: GitBranch,
      color: 'text-purple-600',
      bgColor: 'bg-purple-50',
    },
    {
      label: 'Pseudotime范围',
      value: results.pseudotime 
        ? `${Math.min(...results.pseudotime).toFixed(3)} - ${Math.max(...results.pseudotime).toFixed(3)}`
        : 'N/A',
      icon: Clock,
      color: 'text-green-600',
      bgColor: 'bg-green-50',
    },
    {
      label: 'Entropy范围',
      value: results.entropy 
        ? `${Math.min(...results.entropy).toFixed(3)} - ${Math.max(...results.entropy).toFixed(3)}`
        : 'N/A',
      icon: Activity,
      color: 'text-blue-600',
      bgColor: 'bg-blue-50',
    },
    {
      label: '分析模式',
      value: results.simulated ? '模拟模式' : 'SAC算法',
      icon: Sparkles,
      color: results.simulated ? 'text-yellow-600' : 'text-green-600',
      bgColor: results.simulated ? 'bg-yellow-50' : 'bg-green-50',
    },
  ]

  return (
    <div className="bg-white rounded-xl shadow-sm border border-gray-100 p-6 animate-fade-in">
      <h2 className="text-lg font-semibold text-gray-800 mb-4 flex items-center gap-2">
        <Sparkles className="w-5 h-5 text-purple-500" />
        分析结果摘要
      </h2>
      
      <div className="grid grid-cols-2 gap-4">
        {stats.map((stat, index) => {
          const Icon = stat.icon
          return (
            <div 
              key={index}
              className={`p-4 rounded-xl ${stat.bgColor} transition-transform hover:scale-105`}
            >
              <div className="flex items-center gap-3">
                <div className={`p-2 rounded-lg bg-white shadow-sm`}>
                  <Icon className={`w-5 h-5 ${stat.color}`} />
                </div>
                <div>
                  <p className="text-sm text-gray-600">{stat.label}</p>
                  <p className={`text-lg font-semibold ${stat.color}`}>{stat.value}</p>
                </div>
              </div>
            </div>
          )
        })}
      </div>

      {results.simulated && (
        <div className="mt-4 p-3 bg-yellow-50 border border-yellow-200 rounded-lg">
          <p className="text-sm text-yellow-700">
            ⚠️ 当前使用模拟数据模式（未安装scFocus库）。实际分析需要在后端环境中安装scFocus。
          </p>
        </div>
      )}
    </div>
  )
}
