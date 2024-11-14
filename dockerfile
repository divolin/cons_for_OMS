FROM python:3.9 AS builder

WORKDIR /app

COPY . /app

RUN pip install --no-cache-dir biopython==1.83 cutadapt==4.8 Flask

RUN chmod +x /app/muscle3.8.31_i86linux64
RUN chmod +x /app/minimap2-2.28_x64-linux/minimap2

FROM python:3.9-slim

WORKDIR /app
COPY --from=builder /app /app
COPY --from=builder /usr/local/lib/python3.9/site-packages /usr/local/lib/python3.9/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin 

ENV PATH="/app:${PATH}"
ENV FLASK_APP=app
ENV FLASK_RUN_HOST=0.0.0.0

EXPOSE 5000

CMD ["python", "app.py"]
