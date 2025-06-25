#!/bin/bash

# Usage: ./script.sh dir1 dir2 [dir3 ...]

if [ "$#" -lt 2 ]; then
  echo "You must provide paths to at least two directories."
  exit 1
fi

# make outdir if it doesn't exist
output_dir="in_fq"
mkdir -p "$output_dir"

# define subscripts
index_to_subscript() {
    local index=$1
    local sub=""
    local char

    while true; do
      char=$((index % 26))
      sub=$(printf "\\$(printf '%03o' $((char + 65)))")$sub
      index=$((index / 26 - 1))
      [[ $index -lt 0 ]] && break
    done

    echo "$sub"
}

# verify dirs
for i in "$@"; do
  if [ ! -d "$i" ]; then
    echo "Directory $i does not exist."
    exit 1
  fi
done

# process each dir
index=0
for dir in "$@"; do
  subscript=$(index_to_subscript $index)

  for file in "$dir"/*; do
    if [ -f "$file" ]; then
      target=$(readlink -f "$file")
      filename=$(basename "$file")

      name="${filename%%.*}"
      extension="${filename#${name}}"

      symlink_name="${output_dir}/${name}_${subscript}${extension}"

      ln -s "$target" "$symlink_name"
    fi
  done

  index=$((index + 1))
done

echo "Symlinks created in the '$output_dir' directory."
