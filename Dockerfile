FROM ubuntu:20.04
RUN apt-get update && apt-get -y update
RUN apt-get install -y build-essential python3.8 python3-pip python3-dev
RUN pip3 -q install pip --upgrade
RUN apt-get install -y python3-venv



RUN mkdir exchange
RUN mkdir clusty


ADD ./clusty/ /clusty

ENV MODULE2_4=/venv/module2_4
RUN python3 -m venv $MODULE2_4
ENV PATH="$MODULE2_4/bin:$PATH"

COPY module2_4requirements.txt .
RUN pip install -r module2_4requirements.txt

ENV MODULE1_3=/venv/module1_3
RUN python3 -m venv $MODULE1_3
ENV PATH="$MODULE1_3/bin:$PATH"

COPY module1_3requirements.txt .
RUN pip install -r module1_3requirements.txt
RUN rm module2_4requirements.txt
RUN rm module1_3requirements.txt
