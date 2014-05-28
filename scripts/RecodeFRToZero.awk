 {
	if ($4 == "-") {
		print $1,$2,$3,"0","0"
	} else {
		print $0
	}
}