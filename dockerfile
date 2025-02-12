FROM pypy:3.9-slim AS builder

WORKDIR /app

RUN apt-get update && apt-get install -y \
    build-essential \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libffi-dev \
    git \

    && apt-get clean

COPY . /app


RUN pip install --no-cache-dir biopython==1.83 cutadapt==4.8 Flask

RUN chmod +x /app/muscle3.8.31_i86linux64
RUN chmod +x /app/minimap2-2.28_x64-linux/minimap2


FROM pypy:3.9-slim

WORKDIR /app


RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    && apt-get clean

COPY --from=builder /app /app
COPY --from=builder /opt /opt

ENV PATH="/app:${PATH}"
ENV FLASK_APP=app
ENV FLASK_RUN_HOST=0.0.0.0

EXPOSE 5000

CMD ["sh", "-c", "umask 000 && pypy app.py"]

