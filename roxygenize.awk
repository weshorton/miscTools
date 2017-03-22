#!/bin/awk -f

BEGIN{
    in_function = "false";
}
{
#print($0)
    if ($0~/^function/)
    {
	        #print "has function";
	in_function = "true";
    }

    if (in_function == "true")
    {
	        #print "in function"
	gsub(/# /, "#\047 ");
    }

    if (NF == 0)
    {
	        #print "blank line";
	in_function = "false";
    }

        #print in_function
        print $0
}

BEGIN{
    in_function = "false";
}
{
    # print($0)
    if ($0~/^function/)
    {
	# print "has function";
	in_function = "true";
    }

    if (in_function == "true");
    {
	# print "in function"
	gsub(/# /, "#\047 ");
	gsub(/# [a-zA-z]+
