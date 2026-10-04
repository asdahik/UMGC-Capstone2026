FROM python:3.13-slim

ENV CUDA_VISIBLE_DEVICES=""

RUN pip install --no-cache-dir torch==2.14.0 --index-url https://download.pytorch.org/whl/cpu && \
    pip install --no-cache-dir \
      casanovo==5.2.1 \
      pyteomics==5.0.1 \
      lxml \
      pandas

WORKDIR /app
CMD ["/bin/bash"]