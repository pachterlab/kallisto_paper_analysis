#!/bin/bash

rsync -rv --include '*/' --include '*.h5' --include '*.sf' --include '*isoforms.results' --include 'results.xprs' --exclude '*' --prune-empty-dirs . 30000000_summary
