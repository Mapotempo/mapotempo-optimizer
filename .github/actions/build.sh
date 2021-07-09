#!/usr/bin/env bash

set -e

if [[ $GITHUB_REF == *tags* ]]; then
  TAG=${GITHUB_REF#refs/tags/}
else
  echo "> Building image for dev purpose (${REPOSITORY}) with or-tools ${ORTOOLS_VERSION}"
  TAG=$(echo "${REPOSITORY}" | sed "s/\/.*//")
fi

IMAGE_NAME=${REGISTRY}/mapotempo-ce/optimizer-ortools:${TAG}

case $ORTOOLS_VERSION in
  'v7.1')
    export ORTOOLS_URL="https://github.com/google/or-tools/releases/download/v7.1/or-tools_debian-9_v7.1.6720.tar.gz"
    ;;
  'v7.4')
    export ORTOOLS_URL="https://github.com/google/or-tools/releases/download/v7.4/or-tools_debian-9_v7.4.7247.tar.gz"
    ;;
  'v7.5')
    export ORTOOLS_URL="https://github.com/google/or-tools/releases/download/v7.5/or-tools_debian-10_v7.5.7466.tar.gz"
    ;;
  'v7.8')
   export ORTOOLS_URL="https://github.com/google/or-tools/releases/download/v7.8/or-tools_debian-10_v7.8.7959.tar.gz"
    ;;
  *)
    echo "Unknown OR-Tools version"
esac

echo "Download asset at ${ORTOOLS_URL}"
docker build --build-arg ORTOOLS_URL=${ORTOOLS_URL} -f ./Dockerfile -t "${IMAGE_NAME}" .
docker run -d --name optimizer -t "${IMAGE_NAME}"
docker exec -i optimizer bash -c "LD_LIBRARY_PATH=/srv/or-tools/lib/ /srv/optimizer-ortools/tsp_simple -time_limit_in_ms 500 -intermediate_solutions -instance_file '/srv/optimizer-ortools/data/49Missions_7Vehicles_VRP2TW' -solution_file '/tmp/optimize-or-tools-output'"

echo "::set-output name=image_name::${IMAGE_NAME}"
