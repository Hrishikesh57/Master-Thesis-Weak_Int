FROM docker.io/library/python:3.14-slim-bookworm
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1 MPLCONFIGDIR=/tmp/matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends gfortran gcc binutils liblapack-dev libblas-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements-lock.txt .
RUN python -m pip install --no-cache-dir -r requirements-lock.txt
COPY workspace ./workspace
COPY scripts ./scripts
RUN python workspace/programs/weak_spherical/build.py && python -m unittest discover -s workspace/programs/weak_spherical -p 'test_spherical.py' -v
RUN mkdir /output && chmod 1777 /output
USER 65532:65532
ENTRYPOINT ["python", "scripts/reproduce.py"]
CMD ["--output", "/output/demo"]
