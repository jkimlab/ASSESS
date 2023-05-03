#!/bin/bash
echo Docker build start
date
docker build -t assessw . > log.dockerbuild.txt 2>&1
date
echo Finished
