'use client'

import { createContext, useContext, useState, useEffect, ReactNode } from 'react'
import { translations, Language } from '@/lib/i18n'

interface LanguageContextType {
  lang: Language
  setLang: (lang: Language) => void
  t: typeof translations.en
}

const LanguageContext = createContext<LanguageContextType | undefined>(undefined)

export function LanguageProvider({ children }: { children: ReactNode }) {
  const [lang, setLangState] = useState<Language>('zh')

  useEffect(() => {
    // 从 localStorage 读取语言设置
    const savedLang = localStorage.getItem('scfocus-lang') as Language
    if (savedLang && (savedLang === 'en' || savedLang === 'zh')) {
      setLangState(savedLang)
    } else {
      // 根据浏览器语言自动设置
      const browserLang = navigator.language.toLowerCase()
      setLangState(browserLang.startsWith('zh') ? 'zh' : 'en')
    }
  }, [])

  const setLang = (newLang: Language) => {
    setLangState(newLang)
    localStorage.setItem('scfocus-lang', newLang)
  }

  const t = translations[lang]

  return (
    <LanguageContext.Provider value={{ lang, setLang, t }}>
      {children}
    </LanguageContext.Provider>
  )
}

export function useLanguage() {
  const context = useContext(LanguageContext)
  if (!context) {
    throw new Error('useLanguage must be used within a LanguageProvider')
  }
  return context
}
