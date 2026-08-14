FROM docker.io/library/python:3.12-bookworm

COPY . .

RUN echo ls

COPY ./requirements.txt /src/requirements.txt

WORKDIR /src

# RUN pip install -r requirements.txt 

COPY src/* /src/

RUN pip install .

CMD [ "python", "main.py"]