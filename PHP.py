def predictVirusHost(scriptPath,bacteriaKmerDir,bacteriaKmerName,outFileDir,dicVirusSeqLength):
    import pandas as pd
    import numpy  as np
    import joblib
    
    modelFullLength = joblib.load(scriptPath+'/model/FullLength/FullLength.m')
    model10k = joblib.load(scriptPath+'/model/10k/10k.m')
    model5k = joblib.load(scriptPath+'/model/5k/5k.m')
    model3k = joblib.load(scriptPath+'/model/3k/3k.m')
    model1k = joblib.load(scriptPath+'/model/1k/1k.m')
    ###
    hostAll = pd.read_csv(bacteriaKmerDir + bacteriaKmerName, sep=',', header=None, index_col=0).astype('float32')  # 全部的细菌基因组kmer
    listHostName = hostAll._stat_axis.values.tolist()  # 得到测试的宿主的名字列表
    print('bacteriaNum', len(listHostName))
    fileOut = open(outFileDir + bacteriaKmerName + '_Prediction_Maxhost.csv', 'w')
    fileOut.write( 'queryVirus\tscore\t_maxScoreHost\n')
    fileOutAll = open(outFileDir + bacteriaKmerName + '_Prediction_Allhost.csv','w')
    fileOutAll.write('host')
    for eachHost in listHostName:
        fileOutAll.write(','+str(eachHost))
    fileOutAll.write('\n')
    
    testVirusAll = pd.read_csv(outFileDir+'virusKmer', sep=',', header=None, index_col=0)
    testList = testVirusAll._stat_axis.values.tolist()
    n = 0
    for eachVirus in testList:
        n+=1
        print('Counting score\t',eachVirus,str(n)+'/'+str(len(testList)))
        length = int(dicVirusSeqLength[eachVirus])
        if length >12500:
            model = modelFullLength
        elif length>7500 and length<=12500:
            model = model10k
        elif length>4000 and length<=7500:
            model = model5k
        elif length>2000 and length<=4000:
            model = model3k
        elif length>4 and length<=2000:
            model = model1k
        else:
            continue
        
        virusParameter = testVirusAll.loc[eachVirus]
        temp1 = hostAll.sub((virusParameter), axis=1)
        temp2 = temp1.sub((temp1), axis=1)
        dataVirusMinusAllHost = temp2.sub((temp1))
        
        pre = model.score_samples(dataVirusMinusAllHost)
        tempMax = 0.0
        for i in range(0,len(pre)):
            if pre[i]>=tempMax: #选取最大的输出
                tempMax = pre[i]
                tempHost= listHostName[i]
        
        fileOut.write( str(eachVirus) +'\t' + str(tempMax) + '\t' + str(tempHost)  + '\n')
        
        #outputAll
        tempMax = 0.0
        fileOutAll.write(eachVirus)
        for i in range(0,len(pre)):
            fileOutAll.write(','+str(pre[i]))
        fileOutAll.write('\n')
    fileOut.close()
    fileOutAll.close()


def main():
    import pandas as pd
    import numpy  as np
    import countKmer
    import datetime,os,sys,getopt

    
    ###default setting
    virusFastaFileDir = ''
    outFileDir = './exampleOutput/'
    bacteriaKmerDir = ''
    bacteriaKmerName = ''
    scriptPath = sys.path[0]+'/'
    print(scriptPath)
    
    opts, args = getopt.getopt(sys.argv[1:], "hv:h:o:c:",["virusFastaFileDir=","outFileDir=","bacteriaKmerDir=","bacteriaKmerName="])
    for op, value in opts:
        if   op == "--virusFastaFileDir":
            virusFastaFileDir = value+'/'
        elif op == "--outFileDir":
            outFileDir = value+'/'
        elif op == "--bacteriaKmerDir":
            bacteriaKmerDir = value+'/'
        elif op == "--bacteriaKmerName":
            bacteriaKmerName = value


    if not os.path.exists(outFileDir):
        os.mkdir(outFileDir)
    if (not os.path.exists(virusFastaFileDir)) or virusFastaFileDir == '':
        print('please input a correct virusFastaFileDir.\ndone.')
        return
    if (not os.path.exists(bacteriaKmerDir)) or bacteriaKmerDir == '':
        print('please input a correct bacteriaKmerDir.\ndone.')
        return
    if (not os.path.exists(bacteriaKmerDir+bacteriaKmerName)) or bacteriaKmerName == '':
        print(bacteriaKmerName+' is not exist in '+bacteriaKmerDir+'.\ndone.')
        return


    dicVirusSeqLength = countKmer.getKmer(virusFastaFileDir,outFileDir,'virusKmer')
    predictVirusHost(scriptPath,bacteriaKmerDir,bacteriaKmerName,outFileDir,dicVirusSeqLength)
    

if __name__ == '__main__':
    main()
    #python3 countKmer.py --fastaFileDir  ./exampleHostGenome --kmerFileDir ./exampleOutput --kmerName HostKmer
    #python3 PHP.py --virusFastaFileDir ./exampleVirusGenome  --outFileDir ./exampleOutput  --bacteriaKmerDir ./exampleOutput  --bacteriaKmerName HostKmer
    
    
    
    
    
    
    
    
    

