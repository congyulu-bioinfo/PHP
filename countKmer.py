def getKmer(dirVirus,outDir,outName):
    import os
    fileOut = open(outDir+'/'+outName,'w')

    virusNameList = os.listdir(dirVirus)
    N = len(virusNameList)
    dicLength = {}
    for n in range(N):
        print('Counting kmer\t'+str(n+1)+'/'+str(N))
        eachVirusName = virusNameList[n]
        fileIn = open(dirVirus +'/'+ eachVirusName,'r')
        textIn = fileIn.readlines()
        fileIn.close()
        seq = ''
        for eachLine in textIn:
            if '>' not in eachLine:
                seq += eachLine.strip('\n')
        dicts, seqLength = countKmerNum(seq)
        dicLength[eachVirusName] = seqLength
        
        seqLength = float(seqLength)
        fileOut.write(eachVirusName)
        for eachItem in sorted(dicts.items(), key=lambda e: e[0]):
            eachValue = eachItem[1]/(seqLength - 1)
            fileOut.write(',' + '%.6f' % eachValue)
        fileOut.write('\n')
        
    fileOut.close()
    return dicLength
    
def countKmerNum(seq):
    seq = seq.upper()
    dic_kmer4 = eval("{'AAGC': 0, 'TTTC': 0, 'ACAG': 0, 'CAGA': 0, 'CAGG': 0, 'GAGA': 0, 'AACT': 0, 'ATAT': 0, 'GCTA': 0, 'CTCA': 0, 'CAGT': 0, 'TGAA': 0, 'AGTC': 0, 'ACAA': 0, 'AAGG': 0, 'CGTG': 0, 'TGAT': 0, 'CTTA': 0, 'CGAC': 0, 'AGGG': 0, 'CACT': 0, 'TGGA': 0, 'CGCA': 0, 'TGGC': 0, 'CCTT': 0, 'TACA': 0, 'AGCT': 0, 'GACT': 0, 'AAAT': 0, 'TAGG': 0, 'TACT': 0, 'TCAG': 0, 'ATTA': 0, 'CTTC': 0, 'ATGG': 0, 'CTGG': 0, 'TCTT': 0, 'TTAT': 0, 'ATGT': 0, 'TTCT': 0, 'TTCG': 0, 'CCCC': 0, 'TGAC': 0, 'AATC': 0, 'TTCA': 0, 'AGGA': 0, 'TTCC': 0, 'CGCG': 0, 'GGTC': 0, 'CAGC': 0, 'GTCC': 0, 'AAAC': 0, 'GCCA': 0, 'TCCG': 0, 'TAGC': 0, 'ATCA': 0, 'ACTG': 0, 'CACG': 0, 'ACTC': 0, 'ATAG': 0, 'ACAT': 0, 'GACG': 0, 'TCTG': 0, 'TCGG': 0, 'GCCT': 0, 'CAAA': 0, 'CCGT': 0, 'CGAA': 0, 'GAAG': 0, 'GTCT': 0, 'CGGG': 0, 'TATG': 0, 'CAAT': 0, 'TAAT': 0, 'CTTT': 0, 'GGAC': 0, 'TGGG': 0, 'GCAT': 0, 'AGCA': 0, 'TCTA': 0, 'GTCA': 0, 'CATG': 0, 'TAAC': 0, 'TAGT': 0, 'TAAG': 0, 'TTTT': 0, 'CCCT': 0, 'TCCA': 0, 'TAAA': 0, 'AGAA': 0, 'TGCG': 0, 'GTAG': 0, 'CACC': 0, 'TCAC': 0, 'AAAA': 0, 'AAGT': 0, 'AACA': 0, 'GTCG': 0, 'TCCT': 0, 'ACTA': 0, 'CTAC': 0, 'GATC': 0, 'CACA': 0, 'ACCG': 0, 'GTAC': 0, 'GTTC': 0, 'TAGA': 0, 'TGCA': 0, 'AGCC': 0, 'TTTA': 0, 'GCAC': 0, 'ATGA': 0, 'AACC': 0, 'CTGT': 0, 'ACTT': 0, 'CTAT': 0, 'AAGA': 0, 'GCGC': 0, 'CGTC': 0, 'CCAG': 0, 'TGTC': 0, 'AATT': 0, 'CGTA': 0, 'GTTT': 0, 'CGGA': 0, 'TCGA': 0, 'TTGT': 0, 'GTAA': 0, 'CCGC': 0, 'GCTT': 0, 'ATCG': 0, 'GAAC': 0, 'GCGA': 0, 'ATTC': 0, 'CGAG': 0, 'TCCC': 0, 'GATA': 0, 'TTGA': 0, 'GCAA': 0, 'AGTG': 0, 'TCTC': 0, 'TGTG': 0, 'ATGC': 0, 'GGTA': 0, 'TGAG': 0, 'GCGG': 0, 'ACCT': 0, 'GAAT': 0, 'CTGA': 0, 'GCGT': 0, 'AGAC': 0, 'GCTG': 0, 'GTGT': 0, 'ATCT': 0, 'CGGC': 0, 'CCCA': 0, 'TGCC': 0, 'CTGC': 0, 'AGCG': 0, 'CCAA': 0, 'CATT': 0, 'CTCT': 0, 'GGTT': 0, 'TCGT': 0, 'GGGC': 0, 'ATAA': 0, 'CGAT': 0, 'TATT': 0, 'ACGC': 0, 'CGCT': 0, 'TATC': 0, 'TCGC': 0, 'CTCC': 0, 'CTAG': 0, 'GGAA': 0, 'AGAT': 0, 'ATAC': 0, 'CCTA': 0, 'CGTT': 0, 'GGTG': 0, 'CCGA': 0, 'AAAG': 0, 'CCGG': 0, 'GGCT': 0, 'ATTT': 0, 'GTGA': 0, 'GGCA': 0, 'TTAG': 0, 'GGGA': 0, 'GCCC': 0, 'GTGG': 0, 'GCAG': 0, 'GTTG': 0, 'GAAA': 0, 'GTTA': 0, 'CGCC': 0, 'TTGG': 0, 'GAGC': 0, 'CTCG': 0, 'AGGT': 0, 'TACC': 0, 'ATTG': 0, 'AATG': 0, 'CAAG': 0, 'AGTT': 0, 'ACCC': 0, 'CCAC': 0, 'CTTG': 0, 'TTAC': 0, 'GAGG': 0, 'GGCG': 0, 'TCAA': 0, 'AGAG': 0, 'CAAC': 0, 'CCTC': 0, 'GTAT': 0, 'ACGA': 0, 'ACGT': 0, 'AATA': 0, 'TCAT': 0, 'GGAT': 0, 'GGAG': 0, 'TGTT': 0, 'AACG': 0, 'ACGG': 0, 'GATT': 0, 'GACC': 0, 'AGGC': 0, 'CGGT': 0, 'CCTG': 0, 'TGTA': 0, 'GAGT': 0, 'GGCC': 0, 'GCCG': 0, 'CATC': 0, 'ACAC': 0, 'GTGC': 0, 'TGCT': 0, 'GGGG': 0, 'ACCA': 0, 'TTAA': 0, 'AGTA': 0, 'GACA': 0, 'TGGT': 0, 'CTAA': 0, 'GGGT': 0, 'TACG': 0, 'TATA': 0, 'CCCG': 0, 'CATA': 0, 'TTTG': 0, 'TTGC': 0, 'GCTC': 0, 'CCAT': 0, 'ATCC': 0, 'GATG': 0}")

    seqLength = len(seq)
    for locus in range(0,seqLength):
        if locus + 4 <= seqLength:
            kmer4 = seq[locus:locus+4]
            if kmer4 in dic_kmer4: 
                dic_kmer4[kmer4] += 1
    return dic_kmer4,seqLength



if __name__ == '__main__':
    import datetime
    import os
    import datetime,os,sys,getopt
    
    outFileDir = './exampleOutput/'
    outName = 'hostkmer'
    bacteriaFastaFileDir = ''
    bacteriaKmerName = ''
    
    ###user settings
    opts, args = getopt.getopt(sys.argv[1:], "hf:d:n:",["help","fastaFileDir=","kmerFileDir=","kmerName="])
    for op, value in opts:
        if op == "--fastaFileDir" or op == "-f":
            bacteriaFastaFileDir = value+'/'
        elif op == "--kmerFileDir" or op == "-d":
            outFileDir = value+'/'
        elif op == "--kmerName" or op == "-n":
            bacteriaKmerName = value
        elif op == "--help" or op == "-h":
            print('Step 1: calculate the K-mer frequency of the host\n')
            print('    python3 countKmer.py --fastaFileDir  ./exampleHostGenome --kmerFileDir ./exampleOutput --kmerName HostKmer\n')
            print('Or use the simplify command\n')
            print('    python3 countKmer.py -f ./exampleHostGenome -d ./exampleOutput -n HostKmer\n')
            print('--fastaFileDir or -f: The fasta file of prokaryotic genome sequences, one genome per file.')
            print('--kmerFileDir or -d: The path of prokaryotic K-mer file.')
            print('--kmerName or -n: The name of prokaryotic K-mer file.\n')
            sys.exit()
            
    if not os.path.exists(outFileDir):
        os.mkdir(outFileDir)
    if bacteriaKmerName == '':
        print('please input a bacteriaKmerName.\nbreak')
    if (not os.path.exists(bacteriaFastaFileDir)) or bacteriaFastaFileDir == '':
        print('please input a correct bacteriaFastaFileDir.\nbreak')
    else:
        print('counting kmer ...')
        
        getKmer(bacteriaFastaFileDir,outFileDir,bacteriaKmerName)
    print('done.')
    
    #python3 countKmer.py --fastaFileDir  ./exampleHostGenome --kmerFileDir ./exampleOutput --kmerName HostKmer
    #python3 countKmer.py -f ./exampleHostGenome -d ./exampleOutput -n HostKmer
    
    
    
    
    
    
    
    
    
