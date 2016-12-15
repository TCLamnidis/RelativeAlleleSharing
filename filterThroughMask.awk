BEGIN {
    pos_field=2
}
{
    pos = $pos_field
    line = $0
	if($1=="#CHROM") print $0
    else {
	while(mask_end <= pos) {
        if( (getline < maskFile) > 0) {
            mask_begin = $2
            mask_end = $3
        }
        else
            exit 0
    }
    if(pos >= mask_begin && pos < mask_end)
        print line
    if(pos >= mask_end)
        exit 0 
	}
}
END {
    exit 0
}