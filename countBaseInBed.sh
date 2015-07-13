args=("$@")
awk '{total_bases+=$3-$2} END {printf "%s\n", total_bases}' ${args[0]}
