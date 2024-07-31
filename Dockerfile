FROM python:3.11

COPY . .

RUN python -m pip install --upgrade pip \
    && pip install . \
    && mclumi --help .

WORKDIR .