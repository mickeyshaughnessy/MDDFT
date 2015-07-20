#!/bin/bash
sqlite3 $1 < check_db.sqlite  | tee $1.text
