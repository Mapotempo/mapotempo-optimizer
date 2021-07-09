# 1.0.13 is the latest version containing bundler 2 required for optimizer-api
FROM phusion/passenger-ruby25:1.0.13

ARG ORTOOLS_URL=${ORTOOLS_URL}

LABEL maintainer="Mapotempo <tech@mapotempo.com>"

WORKDIR /srv/

RUN apt-get update > /dev/null && \
  apt-get -y install git wget pkg-config build-essential cmake autoconf libtool zlib1g-dev lsb-release > /dev/null

ADD . /srv/or-tools

RUN wget -qO- $ORTOOLS_URL | tar xz --strip-components=1 -C /srv/or-tools

ADD . /srv/optimizer-ortools

WORKDIR /srv/optimizer-ortools
RUN make tsp_simple
