#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
image_path="${SOMATIX_SIF:-${repo_dir}/somatix.sif}"
runtime="${SINGULARITY_CMD:-}"

if [[ -z "${runtime}" ]]; then
    if command -v apptainer >/dev/null 2>&1; then
        runtime="apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        runtime="singularity"
    else
        echo "ERROR: neither apptainer nor singularity was found in PATH" >&2
        exit 1
    fi
fi

if [[ ! -f "${image_path}" ]]; then
    echo "ERROR: missing SomatiX image: ${image_path}" >&2
    echo "Build it with: ${script_dir}/build_somatix_sif.sh ${image_path}" >&2
    exit 1
fi

exec "${runtime}" run \
    --bind "${PWD}:${PWD}" \
    --pwd "${PWD}" \
    "${image_path}" \
    "$@"
