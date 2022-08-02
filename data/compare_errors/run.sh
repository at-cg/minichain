#!/bin/bash
sh extract_wrong_mappings.sh Linear.gaf
sh extract_wrong_mappings.sh 10H.gaf
sh extract_wrong_mappings.sh 95H.gaf
python3 compare_erros.py
