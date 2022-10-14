# Replace space by _ in file names
find  -name *.xlsx -type f -print0 | while read -d $'\0' f; do mv -v  ; done
# Rename, organise files by exp type and convert to csv 
./rename_files.sh ./*/*/*/*.xlsx
