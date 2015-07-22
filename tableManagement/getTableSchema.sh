#!/bin/bash

if (( "$#" < 1 )); then
    echo "missing argument: Dataset.TableName"
fi

bq --format=prettyjson show $1 > table.json
sed -n '/fields/,$p' table.json | sed 's/"fields": //' | sed '/]/q' > newTable.json

echo "You might want to use:
    bq update $1 newTable.json
to update the table once you have modified it
"
