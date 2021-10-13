#!/usr/bin/env bash

GIT_DIR=$(git rev-parse --git-dir)

if [ ! -L ${GIT_DIR}/hooks/pre-commit ]; then
  echo "Installing hooks..."
  ln -s ../../scripts/pre-commit.bash $GIT_DIR/hooks/pre-commit
  echo "Done!"
fi
