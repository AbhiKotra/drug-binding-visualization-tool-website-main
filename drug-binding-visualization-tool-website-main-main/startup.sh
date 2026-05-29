#!/bin/bash
gunicorn app:app --timeout 600 --workers 1 --threads 4
