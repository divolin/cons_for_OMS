
FROM python:3.9


WORKDIR /app


RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*


RUN git clone --branch flask_main https://github.com/divolin/cons_for_OMS.git .


RUN pip install --no-cache-dir biopython==1.83 cutadapt==4.8 Flask Werkzeug


RUN chmod +x muscle3.8.31_i86linux64
RUN chmod +x minimap2-2.28_x64-linux/minimap2


ENV PATH="/app:${PATH}"

ENV FLASK_APP=app
ENV FLASK_RUN_HOST=0.0.0.0


EXPOSE 5000


CMD ["python", "app.py"]
