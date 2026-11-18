#!/usr/bin/env sh
# Build the SourceEstimation.jl documentation into docs/_site.
#
#   ./build.sh            build the site
#   ./build.sh --force    re-execute all the Julia chunks
#   ./build.sh --serve    build, then preview with live reload
set -e
cd "$(dirname "$0")"

uv sync --quiet                                                  # Python bridge to the IJulia kernel
julia --project=. -e 'using Pkg; Pkg.instantiate()' >/dev/null   # Julia environment of the docs
julia --project=. _lib/make_api.jl                               # API reference from the docstrings

if [ "$1" = "--serve" ]; then
    calepin watch . _site --serve --open
else
    calepin compile . "$@"
fi
