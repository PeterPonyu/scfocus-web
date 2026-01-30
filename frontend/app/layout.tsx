import type { Metadata } from 'next'
import './globals.css'
import { Providers } from './providers'

export const metadata: Metadata = {
  title: 'scFocus Web - 单细胞强化学习聚焦分析平台',
  description: '使用Soft Actor-Critic算法进行单细胞数据谱系分支识别的Web应用',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="zh-CN">
      <body className="font-sans antialiased">
        <Providers>
          <div className="min-h-screen bg-gradient-to-br from-green-50 via-white to-purple-50">
            {children}
          </div>
        </Providers>
      </body>
    </html>
  )
}
