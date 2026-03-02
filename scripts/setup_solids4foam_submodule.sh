#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."

if [ ! -d external/solids4foam/.git ]; then
    git submodule add https://github.com/solids4foam/solids4foam.git external/solids4foam
fi

git -C external/solids4foam fetch --tags
git -C external/solids4foam checkout v2.3

echo "solids4foam submodule is pinned to $(git -C external/solids4foam describe --tags --always)"
