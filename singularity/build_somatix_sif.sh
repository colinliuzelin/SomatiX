#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"

image_arg="${1:-${repo_dir}/somatix.sif}"
image_dir="$(dirname "${image_arg}")"
image_base="$(basename "${image_arg}")"
image_path="$(cd "${image_dir}" && pwd)/${image_base}"
definition="${2:-${script_dir}/SomatiX.def}"
runtime="${SINGULARITY_CMD:-}"
build_tmp_root="${SOMATIX_BUILD_TMPDIR:-/mnt/nvme4t/tmp/somatix_container_build}"
build_tmp="${build_tmp_root}/tmp"
build_cache="${build_tmp_root}/cache"

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

mkdir -p "$(dirname "${image_path}")"
mkdir -p "${build_tmp}" "${build_cache}"

export TMPDIR="${build_tmp}"
export APPTAINER_TMPDIR="${build_tmp}"
export APPTAINER_CACHEDIR="${build_cache}"
export SINGULARITY_TMPDIR="${build_tmp}"
export SINGULARITY_CACHEDIR="${build_cache}"

echo "[INFO] Repository: ${repo_dir}"
echo "[INFO] Definition: ${definition}"
echo "[INFO] Output SIF: ${image_path}"
echo "[INFO] Runtime: ${runtime}"
echo "[INFO] Build tmp: ${build_tmp}"
echo "[INFO] Build cache: ${build_cache}"

cd "${repo_dir}"
if [[ "${SOMATIX_BUILD_WITH_SUDO:-0}" == "1" ]]; then
    sudo -E "${runtime}" build --tmpdir "${build_tmp}" ${SINGULARITY_BUILD_FLAGS:-} "${image_path}" "${definition}"
else
    "${runtime}" build --tmpdir "${build_tmp}" ${SINGULARITY_BUILD_FLAGS:-} "${image_path}" "${definition}"
fi

echo "[DONE] Built ${image_path}"
echo "[TEST] ${runtime} run ${image_path} --version"
"${runtime}" run "${image_path}" --version
