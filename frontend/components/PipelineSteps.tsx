'use client'

import { Check } from 'lucide-react'
import { LucideIcon } from 'lucide-react'

interface Step {
  id: string
  label: string
  icon: LucideIcon
}

interface PipelineStepsProps {
  steps: Step[]
  currentStep: string
}

export default function PipelineSteps({ steps, currentStep }: PipelineStepsProps) {
  const currentIndex = steps.findIndex(s => s.id === currentStep)

  return (
    <div className="mb-8">
      <div className="flex items-center justify-between relative">
        {/* 背景线 */}
        <div className="absolute left-0 right-0 top-1/2 h-1 bg-gray-200 -translate-y-1/2 rounded-full" />
        
        {/* 进度线 */}
        <div 
          className="absolute left-0 top-1/2 h-1 bg-gradient-to-r from-green-500 to-purple-500 -translate-y-1/2 rounded-full transition-all duration-500"
          style={{ width: `${(currentIndex / (steps.length - 1)) * 100}%` }}
        />

        {steps.map((step, index) => {
          const isCompleted = index < currentIndex
          const isCurrent = index === currentIndex
          const Icon = step.icon

          return (
            <div key={step.id} className="relative z-10 flex flex-col items-center">
              <div
                className={`w-12 h-12 rounded-full flex items-center justify-center transition-all duration-300
                  ${isCompleted 
                    ? 'bg-gradient-to-br from-green-500 to-green-600 shadow-lg shadow-green-200' 
                    : isCurrent 
                      ? 'bg-gradient-to-br from-purple-500 to-purple-600 shadow-lg shadow-purple-200 animate-pulse-slow' 
                      : 'bg-white border-2 border-gray-200'
                  }`}
              >
                {isCompleted ? (
                  <Check className="w-6 h-6 text-white" />
                ) : (
                  <Icon className={`w-6 h-6 ${isCurrent ? 'text-white' : 'text-gray-400'}`} />
                )}
              </div>
              <span className={`mt-2 text-sm font-medium ${
                isCompleted || isCurrent ? 'text-gray-800' : 'text-gray-400'
              }`}>
                {step.label}
              </span>
            </div>
          )
        })}
      </div>
    </div>
  )
}
