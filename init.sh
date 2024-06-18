#!/bin/bash

FILENAME="~/.my.cnf"
TEXT="[AAVengeR]
user=aavenger
password=iAmAavengeR1
host=174.129.238.44
port=3306
database=AAVengeR"


if [ ! -e "$FILENAME" ]; then
    touch $FILENAME
    echo "$TEXT" > "$FILENAME"
    echo "Config created and written."
else
    echo "$TEXT" >> "$FILENAME"
    echo "Config added to cnf."
fi
