cat test.sam | cut -f10,14 | awk '
{
    new_methylstring = substr($2,5,length($2))
    start = match(new_methylstring, /[Zz]/)
    end   = match(new_methylstring, /[Zz][^Zz]*$/)
    if (end - start <= 33) end = start + 33 
    new_read = substr($1,start,end)
    # print( new_read )
    print (start ? start : "NaN"), (end ? end : "NaN")
}' > output
