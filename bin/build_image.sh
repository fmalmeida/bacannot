#!/bin/bash
name=$(basename $(pwd))
docker build -t fmalmeida/bacannot:${1}_${name} .
