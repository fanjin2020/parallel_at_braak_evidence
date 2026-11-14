#!/usr/bin/env bash
set -euo pipefail
umask 077
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1 PYTHONIOENCODING=utf-8
unset JAD_OFFLINE
exec "${PYTHON:-python3}" -u "$HERE/scripts/download_public.py" "$@"

