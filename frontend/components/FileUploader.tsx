'use client'

import { useCallback, useState } from 'react'
import { Upload, FileUp, X, CheckCircle } from 'lucide-react'
import { useLanguage } from '@/lib/LanguageContext'

interface FileUploaderProps {
  onUpload: (file: File) => void
  isLoading: boolean
}

export default function FileUploader({ onUpload, isLoading }: FileUploaderProps) {
  const { t } = useLanguage()
  const [dragActive, setDragActive] = useState(false)
  const [selectedFile, setSelectedFile] = useState<File | null>(null)

  const handleDrag = useCallback((e: React.DragEvent) => {
    e.preventDefault()
    e.stopPropagation()
    if (e.type === 'dragenter' || e.type === 'dragover') {
      setDragActive(true)
    } else if (e.type === 'dragleave') {
      setDragActive(false)
    }
  }, [])

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault()
    e.stopPropagation()
    setDragActive(false)

    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      const file = e.dataTransfer.files[0]
      setSelectedFile(file)
    }
  }, [])

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files[0]) {
      setSelectedFile(e.target.files[0])
    }
  }

  const handleUpload = () => {
    if (selectedFile) {
      onUpload(selectedFile)
    }
  }

  const clearFile = () => {
    setSelectedFile(null)
  }

  return (
    <div className="space-y-4">
      <div
        className={`relative border-2 border-dashed rounded-xl p-8 text-center transition-all cursor-pointer
          ${dragActive 
            ? 'border-green-500 bg-green-50' 
            : 'border-gray-300 hover:border-green-400 hover:bg-gray-50'
          }`}
        onDragEnter={handleDrag}
        onDragLeave={handleDrag}
        onDragOver={handleDrag}
        onDrop={handleDrop}
        onClick={() => document.getElementById('file-input')?.click()}
      >
        <input
          id="file-input"
          type="file"
          accept=".h5ad,.csv"
          onChange={handleFileChange}
          className="hidden"
        />
        
        <div className="flex flex-col items-center gap-3">
          <div className={`p-3 rounded-full ${dragActive ? 'bg-green-100' : 'bg-gray-100'}`}>
            <Upload className={`w-8 h-8 ${dragActive ? 'text-green-600' : 'text-gray-400'}`} />
          </div>
          <div>
            <p className="text-sm font-medium text-gray-700">
              {t.upload.dragDrop} <span className="text-green-600">{t.upload.clickSelect}</span>
            </p>
            <p className="text-xs text-gray-500 mt-1">
              {t.upload.supportedFormats}
            </p>
          </div>
        </div>
      </div>

      {selectedFile && (
        <div className="flex items-center justify-between p-4 bg-green-50 border border-green-200 rounded-lg animate-fade-in">
          <div className="flex items-center gap-3">
            <CheckCircle className="w-5 h-5 text-green-600" />
            <div>
              <p className="text-sm font-medium text-gray-800">{selectedFile.name}</p>
              <p className="text-xs text-gray-500">
                {(selectedFile.size / (1024 * 1024)).toFixed(2)} MB
              </p>
            </div>
          </div>
          <div className="flex items-center gap-2">
            <button
              onClick={(e) => {
                e.stopPropagation()
                clearFile()
              }}
              className="p-1 text-gray-400 hover:text-gray-600"
            >
              <X className="w-4 h-4" />
            </button>
          </div>
        </div>
      )}

      {selectedFile && (
        <button
          onClick={handleUpload}
          disabled={isLoading}
          className="w-full py-3 bg-gradient-to-r from-green-500 to-green-600 text-white rounded-lg font-medium hover:from-green-600 hover:to-green-700 transition-all disabled:opacity-50 flex items-center justify-center gap-2"
        >
          {isLoading ? (
            <div className="w-5 h-5 border-2 border-white border-t-transparent rounded-full animate-spin" />
          ) : (
            <FileUp className="w-5 h-5" />
          )}
          {isLoading ? t.upload.uploading : t.upload.uploadButton}
        </button>
      )}
    </div>
  )
}
