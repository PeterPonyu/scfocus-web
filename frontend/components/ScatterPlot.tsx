'use client'

import { useEffect, useRef } from 'react'

interface ScatterPlotProps {
  x: number[]
  y: number[]
  colors?: number[]
  colorBy?: string
}

export default function ScatterPlot({ x, y, colors, colorBy }: ScatterPlotProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null)

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    // 设置高DPI
    const dpr = window.devicePixelRatio || 1
    const rect = canvas.getBoundingClientRect()
    canvas.width = rect.width * dpr
    canvas.height = rect.height * dpr
    ctx.scale(dpr, dpr)

    // 清除画布
    ctx.fillStyle = '#fafafa'
    ctx.fillRect(0, 0, rect.width, rect.height)

    // 计算边界
    const padding = 40
    const xMin = Math.min(...x)
    const xMax = Math.max(...x)
    const yMin = Math.min(...y)
    const yMax = Math.max(...y)
    const xRange = xMax - xMin || 1
    const yRange = yMax - yMin || 1

    // 归一化颜色值
    let colorMin = 0, colorMax = 1
    if (colors && colors.length > 0) {
      colorMin = Math.min(...colors)
      colorMax = Math.max(...colors)
    }
    const colorRange = colorMax - colorMin || 1

    // 颜色映射函数 (绿色到紫色渐变)
    const getColor = (value: number) => {
      const normalized = (value - colorMin) / colorRange
      
      // 绿色 -> 黄色 -> 紫色
      const r = Math.round(34 + normalized * (139 - 34))
      const g = Math.round(197 - normalized * (197 - 92))
      const b = Math.round(94 + normalized * (246 - 94))
      
      return `rgba(${r}, ${g}, ${b}, 0.7)`
    }

    // 绘制点
    const width = rect.width - 2 * padding
    const height = rect.height - 2 * padding

    for (let i = 0; i < x.length; i++) {
      const px = padding + ((x[i] - xMin) / xRange) * width
      const py = padding + height - ((y[i] - yMin) / yRange) * height
      
      const color = colors && colors[i] !== undefined 
        ? getColor(colors[i]) 
        : 'rgba(34, 197, 94, 0.6)'

      ctx.beginPath()
      ctx.arc(px, py, 2.5, 0, Math.PI * 2)
      ctx.fillStyle = color
      ctx.fill()
    }

    // 绘制坐标轴
    ctx.strokeStyle = '#e5e7eb'
    ctx.lineWidth = 1
    
    // X轴
    ctx.beginPath()
    ctx.moveTo(padding, rect.height - padding)
    ctx.lineTo(rect.width - padding, rect.height - padding)
    ctx.stroke()
    
    // Y轴
    ctx.beginPath()
    ctx.moveTo(padding, padding)
    ctx.lineTo(padding, rect.height - padding)
    ctx.stroke()

    // 标签
    ctx.fillStyle = '#6b7280'
    ctx.font = '12px Inter, sans-serif'
    ctx.textAlign = 'center'
    ctx.fillText('UMAP1', rect.width / 2, rect.height - 10)
    
    ctx.save()
    ctx.translate(15, rect.height / 2)
    ctx.rotate(-Math.PI / 2)
    ctx.fillText('UMAP2', 0, 0)
    ctx.restore()

    // 颜色图例
    if (colors && colors.length > 0 && colorBy) {
      const legendWidth = 100
      const legendHeight = 12
      const legendX = rect.width - padding - legendWidth
      const legendY = padding

      // 渐变
      const gradient = ctx.createLinearGradient(legendX, 0, legendX + legendWidth, 0)
      gradient.addColorStop(0, 'rgba(34, 197, 94, 0.8)')
      gradient.addColorStop(0.5, 'rgba(234, 179, 8, 0.8)')
      gradient.addColorStop(1, 'rgba(139, 92, 246, 0.8)')
      
      ctx.fillStyle = gradient
      ctx.fillRect(legendX, legendY, legendWidth, legendHeight)
      
      ctx.strokeStyle = '#d1d5db'
      ctx.strokeRect(legendX, legendY, legendWidth, legendHeight)

      // 图例标签
      ctx.fillStyle = '#374151'
      ctx.font = '10px Inter, sans-serif'
      ctx.textAlign = 'left'
      ctx.fillText(colorMin.toFixed(2), legendX, legendY + legendHeight + 12)
      ctx.textAlign = 'right'
      ctx.fillText(colorMax.toFixed(2), legendX + legendWidth, legendY + legendHeight + 12)
      
      ctx.textAlign = 'center'
      ctx.fillText(colorBy.replace('scfocus_', ''), legendX + legendWidth / 2, legendY - 5)
    }

  }, [x, y, colors, colorBy])

  return (
    <div className="w-full aspect-square max-h-[500px] bg-gray-50 rounded-lg overflow-hidden">
      <canvas
        ref={canvasRef}
        className="w-full h-full"
        style={{ display: 'block' }}
      />
    </div>
  )
}
