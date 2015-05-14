import collections

def loadiso():
    f = open('isoforms.txt')
    v = ['isoform']
    for line in f:
        v.append(line.strip())
    return v

def parse_files(fn,iso,cut=0):
    read = None
    val = []
    rsemid = None
    with open(fn) as f:
        for line in f:
            r,tr = line.split()
            if cut > 0:
                r = r[:-cut]
            if read is None:
                read = r
                rsemid = int(read.split('_')[2])
            if r == read:
                if tr != '*':
                    val.append(tr)
            else:
                yield (read,iso[rsemid],val)
                read = r
                rsemid = int(read.split('_')[2])
                val = []
                if tr != '*':
                    val.append(tr)
        
        yield (read,iso[rsemid],val)


def test(fn,cut):
    iso = loadiso()
    x = parse_files(fn,iso,cut)
    numreads=0
    nummapped = 0
    numaccurate=0
    specsum = 0.0
    for r,ids,val in x:
        numreads += 1
        if len(val) > 0:
            nummapped += 1
            if ids in val:
                numaccurate += 1
                specsum += 1.0/len(val)
                
    return numaccurate,nummapped,numreads,specsum

def testboth(fn1,fn2):
    iso = loadiso()
    x = parse_files(fn1,iso,0)
    y = parse_files(fn2,iso,2)

    numr = 0
    nmatch = 0
    nkal = 0
    nbwt = 0
    ndif = 0
    nmissing = 0
    nkal_not_bowtie = 0 

    d = {}
    sensbwt = 0
    for r1,id1,val1 in x:
        d[r1] = set(val1)
        if id1 in val1:
            sensbwt += 1

    

    senskal = 0
    diffkal = 0
    diffbwt = 0

    for r2,id2,val2 in y:
        numr += 1
        s2 = set(val2)
        if id2 in s2:
            senskal += 1

        if r2 in d:
            s1 = d[r2]
            
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
            nkal_not_bowtie += 1

    num_bwt = sum(len(t)>0 for t in d.itervalues())
    nbwt_not_kal = num_bwt - (nmatch+nbwt+nkal+ndif)

    nmissing = len(d) - numr
    numr = len(d)

    return numr, nmatch, nbwt, nkal, ndif, nmissing, nkal_not_bowtie, nbwt_not_kal,diffkal, diffbwt, senskal, sensbwt

    
    

