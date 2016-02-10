import sys
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
                yield (read,iso[rsemid],set(val))
                read = r
                rsemid = int(read.split('_')[2])
                val = []
                if tr != '*':
                    val.append(tr)
        
        yield (read,iso[rsemid],set(val))


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
    y = parse_files(fn2,iso,0)

    r1,r2 = None,None
    id1,id2 = None,None
    s1,s2 = [],[]
    
    numr = 0
    nmatch = 0
    nkal = 0
    nbwt = 0
    ndif = 0
    nmissing = 0
    nkal_not_bowtie = 0 
    nbowtie_not_kal = 0

    sensbwt = 0
    senskal = 0
    diffkal = 0
    diffbwt = 0


    while True:
        xdone,ydone = False, False
	try:
            r1,id1,s1 = x.next()
        except StopIteration:
            xdone=True
        try:
            r2,id2,s2 = y.next()
        except StopIteration:
            ydone = True

        if xdone and ydone:
            break
        if xdone and not ydone:
            print "error, more reads in ",f1,"than",f2
            sys.exit(1)
        if not xdone and ydone:
            print "error, more reads in ",f2,"than",f1
            sys.exit(1)

	numr += 1

        if len(s1)==0 and len(s2)==0:
            nmissing += 1
        elif len(s1)==0 and len(s2) > 0:
            if id2 in s2:
                senskal += 1
            nkal_not_bowtie += 1
        elif len(s1) > 0 and len(s2) ==0:
            if id1 in s1:
                sensbwt += 1
            nbowtie_not_kal += 1
        else:
            if id1 in s1:
                sensbwt +=1
            if id2 in s2:
                senskal += 1
                
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



    print 'Number of read:',numr
    print 'Number of missing:',nmissing
    print 'Number of bowtie missing:',nkal_not_bowtie
    print 'Number of kal missing:',nbowtie_not_kal
    print 'Number of matches:',nmatch
    print 'Number of kal more specific:',nkal
    print 'Number of bowtie more specific:',nbwt
    print 'Bowtie specificity:',sensbwt/float(numr)
    print 'Kallisto specificity:',senskal/float(numr)

    print 'Number of non-simple differences:',ndif
    print 'Average number of kal matches:',(diffkal/float(nbwt+nkal+ndif))
    print 'Average number of bowtie matches:',(diffbwt/float(nbwt+nkal+ndif))
    print 'Fraction of bowtie matches:',(numr-nmissing-nkal_not_bowtie)/float(numr)
    print 'Fraction of kallisto matches:',(numr-nmissing-nbowtie_not_kal)/float(numr)
    print 'Fraction of both matches:',(numr-nmissing-nbowtie_not_kal-nkal_not_bowtie)/float(numr)
    print 'Fraction of exact both matches:',(nmatch)/float(numr-nmissing-nbowtie_not_kal-nkal_not_bowtie)



if __name__ == '__main__':
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    
    testboth(f1,f2)
    

    
    

