'use client'

interface ParamConfig {
  key: string
  label: string
  type: 'number' | 'text' | 'select'
  min?: number
  max?: number
  step?: number
  options?: { value: string | number; label: string }[]
}

interface ParameterPanelProps {
  params: any
  setParams: (params: any) => void
  paramConfig: ParamConfig[]
}

export default function ParameterPanel({ params, setParams, paramConfig }: ParameterPanelProps) {
  const handleChange = (key: string, value: any, type: string) => {
    const newValue = type === 'number' ? parseFloat(value) : value
    setParams({ ...params, [key]: newValue })
  }

  return (
    <div className="space-y-4">
      {paramConfig.map((config) => (
        <div key={config.key} className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            {config.label}
          </label>
          
          {config.type === 'select' ? (
            <select
              value={params[config.key]}
              onChange={(e) => handleChange(config.key, e.target.value, 'text')}
              className="w-full px-3 py-2 bg-gray-50 border border-gray-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-green-500 focus:border-transparent transition-all"
            >
              {config.options?.map((option) => (
                <option key={option.value} value={option.value}>
                  {option.label}
                </option>
              ))}
            </select>
          ) : (
            <input
              type={config.type}
              value={params[config.key]}
              onChange={(e) => handleChange(config.key, e.target.value, config.type)}
              min={config.min}
              max={config.max}
              step={config.step || (config.type === 'number' ? 1 : undefined)}
              className="w-full px-3 py-2 bg-gray-50 border border-gray-200 rounded-lg focus:outline-none focus:ring-2 focus:ring-green-500 focus:border-transparent transition-all"
            />
          )}
          
          {config.min !== undefined && config.max !== undefined && config.type === 'number' && (
            <input
              type="range"
              value={params[config.key]}
              onChange={(e) => handleChange(config.key, e.target.value, config.type)}
              min={config.min}
              max={config.max}
              step={config.step || 1}
              className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-green-500"
            />
          )}
        </div>
      ))}
    </div>
  )
}
