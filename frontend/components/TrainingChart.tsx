'use client'

import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts'

interface TrainingChartProps {
  rewards: number[][] 
  errors: number[][]
}

export default function TrainingChart({ rewards, errors }: TrainingChartProps) {
  // 转换数据格式
  const data = rewards.flatMap((metaRewards, metaIndex) => 
    metaRewards.map((reward, episodeIndex) => ({
      episode: metaIndex * metaRewards.length + episodeIndex,
      reward: reward,
      error: errors[metaIndex]?.[episodeIndex] || 0,
      meta: metaIndex + 1,
    }))
  ).filter((_, i) => i % Math.max(1, Math.floor(rewards.flat().length / 100)) === 0) // 采样以减少点数

  if (data.length === 0) {
    return (
      <div className="h-64 flex items-center justify-center text-gray-500">
        暂无训练数据
      </div>
    )
  }

  return (
    <div className="h-80">
      <ResponsiveContainer width="100%" height="100%">
        <LineChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          <XAxis 
            dataKey="episode" 
            label={{ value: 'Episode', position: 'bottom', offset: -5 }}
            tick={{ fontSize: 12 }}
            stroke="#9ca3af"
          />
          <YAxis 
            yAxisId="left"
            label={{ value: 'Reward', angle: -90, position: 'insideLeft' }}
            tick={{ fontSize: 12 }}
            stroke="#22c55e"
          />
          <YAxis 
            yAxisId="right" 
            orientation="right"
            label={{ value: 'Error', angle: 90, position: 'insideRight' }}
            tick={{ fontSize: 12 }}
            stroke="#8b5cf6"
          />
          <Tooltip 
            contentStyle={{ 
              backgroundColor: 'white', 
              border: '1px solid #e5e7eb',
              borderRadius: '8px',
              boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1)'
            }}
            formatter={(value: number, name: string) => [value.toFixed(4), name]}
          />
          <Legend 
            wrapperStyle={{ paddingTop: '20px' }}
          />
          <Line 
            yAxisId="left"
            type="monotone" 
            dataKey="reward" 
            stroke="#22c55e" 
            strokeWidth={2}
            dot={false}
            name="奖励值"
          />
          <Line 
            yAxisId="right"
            type="monotone" 
            dataKey="error" 
            stroke="#8b5cf6" 
            strokeWidth={2}
            dot={false}
            name="误差值"
          />
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}
