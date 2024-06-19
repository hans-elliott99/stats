
# sudo apt install tree ## if needed
cd ./R
tree -H '.' \
    -L 1 \
    --noreport \
    --dirsfirst \
    --charset utf-8 \
    --ignore-case \
    --timefmt '%d-%b-%Y %H:%M' \
    -I "index.html" \
    -T 'Some Stats Notes' \
    -s -D \
    -P "*.html" \
    -o index.html
