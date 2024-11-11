# Используем базовый образ с Python 3.9
FROM python:3.9

# Устанавливаем рабочую директорию внутри контейнера
WORKDIR /app

# Устанавливаем необходимые системные пакеты
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Клонируем ваш репозиторий с GitHub
RUN git clone https://github.com/divolin/cons_for_OMS.git .

# Устанавливаем необходимые Python-пакеты
RUN pip install --no-cache-dir biopython==1.83 cutadapt==4.8

# Делаем бинарные файлы исполняемыми
RUN chmod +x muscle3.8.31_i86linux64
RUN chmod +x minimap2-2.28_x64-linux/minimap2

# Добавляем директорию с бинарными файлами в PATH
ENV PATH="/app:${PATH}"

# Указываем ENTRYPOINT для запуска скрипта с параметрами
ENTRYPOINT ["python", "consensus.py"]

# CMD оставляем пустым, чтобы можно было передавать аргументы при запуске
CMD []

