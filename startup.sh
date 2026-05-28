#!/bin/bash
pip install -r /home/site/wwwroot/requirements.txt --quiet
gunicorn app:app --timeout 180 --workers 1 --threads 4
