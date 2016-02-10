import sys
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
                yield (read,set(val))
                read = r
                val = []
                if tr != '*':
                    val.append(tr)
        
        yield (read,set(val))



def testboth(fn1,fn2):
    x = parse_files(fn1)
    y = parse_files(fn2)

    numr = 0
    nmatch = 0
    nkal = 0
    nbwt = 0
    ndif = 0
    nkal_not_bowtie = 0 
    nbowtie_not_kal = 0
    n200 = 0
    nmissing = 0
    diffkal = 0
    diffbwt = 0
    r1,r2 = '',''
    s1,s2 = set(),set()

    while True:
        xdone,ydone = False, False
	try:
            r1,s1 = x.next()
        except StopIteration:
            xdone=True
        try:
            r2,s2 = y.next()
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
            nkal_not_bowtie += 1
        elif len(s1)>0 and len(s2)==0:
            if len(s1) > 200:
                n200 += 1
            nbowtie_not_kal += 1
        else:
            if len(s1) > 200:
                n200 += 1
            
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

    nbwt_not_kal = (numr-nmissing-nkal_not_bowtie) - (nmatch+nbwt+nkal+ndif)

    print 'Number of read:',numr
    print 'Number of missing:',nmissing
    print 'Number of bowtie missing:',nkal_not_bowtie
    print 'Number of kal missing:',nbowtie_not_kal
    print 'Number of > 200 matches:',n200
    print 'Number of matches:',nmatch
    print 'Number of kal more specific:',nkal
    print 'Number of bowtie more specific:',nbwt
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
    
