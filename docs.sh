#!/bin/bash
# https://dev.to/deciduously/prepare-your-rust-api-docs-for-github-pages-2n5i
cargo doc --no-deps
rm -rf ./docs
echo "<meta http-equiv=\"refresh\" content=\"0; url=hsc\">" > target/doc/index.html
cp -r target/doc ./docs
