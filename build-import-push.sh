#!/bin/sh

jupyter-book build ../dm4bem_book/ && ghp-import -n -p -f _build/html
