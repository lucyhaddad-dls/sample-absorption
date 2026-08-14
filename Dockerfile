FROM docker.io/library/python:3.12-bookworm

COPY . .

WORKDIR .

RUN pip install .

CMD [ "python", "tst.py"]