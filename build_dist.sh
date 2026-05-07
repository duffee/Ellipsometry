#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

VERSION=$(perl -MExtUtils::MakeMaker -le \
  'print MM->parse_version("lib/Physics/Ellipsometry/VASE.pm")')
DIST="Physics-Ellipsometry-VASE-${VERSION}"

echo "==> Building ${DIST}"

# Clean previous build artifacts
if [[ -f Makefile ]]; then
    make clean >/dev/null 2>&1 || true
fi

# Generate Makefile, build, test, package
perl Makefile.PL
make
make test
make dist

echo ""
echo "==> ${DIST}.tar.gz created successfully"
ls -lh "${DIST}.tar.gz"
