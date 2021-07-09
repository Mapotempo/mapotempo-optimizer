#!/usr/bin/env bash

echo "> Deliver image with or-tools ${ORTOOLS_VERSION}"
echo "$REGISTRY_PASSWORD" | docker login -u "$REGISTRY_USERNAME" --password-stdin "${REGISTRY}"
docker push "${IMAGE_NAME}"
