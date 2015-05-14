import collections


def parse_files(fn):
    read = None
    val = []
    with open(fn) as f:
        for line in f:
            r,tr = line.split()
            if read is None:
                read = r
            if r == read:
                if tr != '*':
                    val.append(tr)
            else:
                yield (read,val)
                read = r
                val = []
                if tr != '*':
                    val.append(tr)
        
        yield (read,val)



def testboth(fn1,fn2):
    x = parse_files(fn1)
    y = parse_files(fn2)

    numr = 0
    nmatch = 0
    nkal = 0
    nbwt = 0
    ndif = 0
    nkal_not_bowtie = 0 
    n200 = 0
    nmissing = 0
    diffkal = 0
    diffbwt = 0

    d = {}
    for r1,val1 in x:
        d[r1] = set(val1)
        if len(val1) >= 200:
            n200 += 1

    for r2,val2 in y:
        s2 = set(val2)

        s1 = d[r2]
        if len(s1) > 0:    
            if s1==s2:
                nmatch += 1
            else:
                if s1 <= s2:
                    nbwt += 1
                elif s2 <= s1:
                    nkal += 1
                else:
                    ndif += 1

                diffbwt += len(s1)
                diffkal += len(s2)
        else:
            if len(s2)>0:
                nkal_not_bowtie += 1
            else:
                nmissing += 1

    num_bwt = sum(len(t)>0 for t in d.itervalues())
    nbwt_not_kal = num_bwt - (nmatch+nbwt+nkal+ndif)

    numr = len(d)

    return numr, nmatch, nbwt, nkal, ndif, nmissing, nkal_not_bowtie, nbwt_not_kal, diffkal, diffbwt, n200

    
    

