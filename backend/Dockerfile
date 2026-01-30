# scFocus-Web 后端 Render 配置
FROM python:3.11-slim

WORKDIR /app

# 安装系统依赖
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 复制 requirements
COPY requirements.txt .

# 安装 Python 依赖
RUN pip install --no-cache-dir -r requirements.txt

# 复制应用代码
COPY . .

# 创建数据目录
RUN mkdir -p /tmp/scfocus_data

# 暴露端口
EXPOSE 10000

# 启动命令
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "10000"]
