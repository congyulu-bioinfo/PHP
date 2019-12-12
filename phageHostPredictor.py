def predictVirusHost(dicVirusSeqLength,testVirusFileName,testVirusFileDir, hostKmerDir, hostKmerName, outFileDir,coreNum):
    import pandas as pd
    import numpy  as np
    from sklearn.externals import joblib
    
    ###fullLengthModel
    coreNum = 1

    dirModelRepeat =  './model/coreNum_'+str(coreNum)+'/'     #modelDir      
    modelFullLength = joblib.load(dirModelRepeat + 'allVirus_randomForest.m')
    ###shortSeqModel
    dicModel = {}
    for each in dicVirusSeqLength:
        length = int(dicVirusSeqLength[each])
        if length>7500 and length<=12500:
            dicModel[10000] = ""
        elif length>4000 and length<=7500:
            dicModel[5000] = ""    
        elif length>2500 and length<=4000:
            dicModel[3000] = "" 
        elif length>1500 and length<=2500:
            dicModel[2000] = ""   
        elif length>750 and length<=1500:
            dicModel[1000] = ""
        elif length>6 and length<=750:
            dicModel[500] = ""
    print(dicModel)
    for length in [500,1000,2000,3000,5000,10000]:
        if 10000 in dicModel and dicModel[10000]=="":
            model10000bp = joblib.load('./model/'+str(10000)+'bp/' + 'd2_randomForest.m')
            dicModel[10000] = "loaded"
        if 5000 in dicModel and dicModel[5000]=="":
            model5000bp = joblib.load('./model/'+str(5000)+'bp/' + 'd2_randomForest.m')  
            dicModel[5000] = "loaded"
        if 3000 in dicModel and dicModel[3000]=="":
            model3000bp = joblib.load('./model/'+str(3000)+'bp/' + 'd2_randomForest.m')  
            dicModel[3000] = "loaded"
        if 2000 in dicModel and dicModel[2000]=="":
            model2000bp = joblib.load('./model/'+str(2000)+'bp/' + 'd2_randomForest.m')    
            dicModel[2000] = "loaded"
        if 1000 in dicModel and dicModel[1000]=="":
            model1000bp = joblib.load('./model/'+str(1000)+'bp/' + 'd2_randomForest.m')    
            dicModel[1000] = "loaded"
        if 500 in dicModel and dicModel[500]=="":
            model500bp = joblib.load('./model/'+str(500)+'bp/' + 'd2_randomForest.m')     
            dicModel[500] = "loaded"
    ###
    hostAll = pd.read_csv(hostKmerDir + hostKmerName, sep=',', header=None, index_col=0).astype('float32')  # 全部的细菌基因组kmer
    print('len(hostAll)', len(hostAll))
    listHostName = hostAll._stat_axis.values.tolist()  # 得到测试的宿主的名字列表
    print(len(listHostName))
    fileOut = open(outFileDir + testVirusFileName + 'Prediction_Maxhost', 'w')
    fileOut.write( 'queryVirus' +'\t' + 'score' + '\t' + 'maxScoreHost\n')
    fileOutAll = open(outFileDir + testVirusFileName + 'Prediction_Allhost','w')
    fileOutAll.write('host')
    for eachHost in listHostName:
        fileOutAll.write('\t'+str(eachHost))
    fileOutAll.write('\n')
    
    testVirusAll = pd.read_csv(testVirusFileDir+testVirusFileName, sep=',', header=None, index_col=0)
    testList = testVirusAll._stat_axis.values.tolist()
    for eachVirus in testList:
        print(eachVirus)
        ###ChooseModelToPredict
        length = int(dicVirusSeqLength[eachVirus])
        print(length)
        if length >12500:
            model = modelFullLength
        elif length>7500 and length<=12500:
            model = model10000bp
            
        elif length>4000 and length<=7500:
            model = model5000bp
                  
        elif length>2500 and length<=4000:
            model = model3000bp
                   
        elif length>1500 and length<=2500:
            model = model2000bp
                 
        elif length>750 and length<=1500:
            model = model1000bp
                 
        elif length>6 and length<=750:
            model = model500bp
        
        virusParameter = testVirusAll.loc[eachVirus]
        temp1 = hostAll.sub((virusParameter), axis=1)
        temp2 = temp1.sub((temp1), axis=1)
        dataVirusMinusAllHost = temp2.sub((temp1))  # 得到该病毒减去所有宿主的数据
        pre = model.predict_proba(dataVirusMinusAllHost)
        #outputMax
        preScoreList = []
        for i in range(0,len(pre)):
            preScoreList.append(pre[i][1])

        tempMax = 0.0
        for i in range(0, len(pre)):
            if preScoreList[i] >= tempMax:  # 选取最大的输出
                tempMax = preScoreList[i]
                tempHost = str(listHostName[i])

        fileOut.write( str(eachVirus) +'\t' + str(tempMax) + '\t' + tempHost  + '\n')
        
        #outputAll
        tempMax = 0.0
        fileOutAll.write(eachVirus)
        for i in range(0,len(pre)):
            fileOutAll.write('\t'+str(pre[i][1]))
        fileOutAll.write('\n')
    fileOut.close()
    fileOutAll.close()

def main():
    import os
    import pandas as pd
    import numpy  as np
    import countKmer
    import datetime,os,sys,getopt

    
    ###default setting
    virusFastaFileDir = ''
    outFileDir = './exampleOutput/'   #outFileDir
    hostFastaFileDir = ''
    coreNum = 1        #computerCoreNumUsed
    dirModelCoreNum =  './model/coreNum_'+str(coreNum)+'/'  #modelDir
    testVirusKmerName = 'virusKmer_1to6_input'             #testVirusKmerName
    ###                                      
    
    opts, args = getopt.getopt(sys.argv[1:], "hv:h:o:c:",["virusFastaFileDir=","hostFastaFileDir=","outFileDir="])
    ###user settings
    for op, value in opts:
        if op == "--virusFastaFileDir":
            virusFastaFileDir = value+'/'
        elif op == "--outFileDir":
            outFileDir = value+'/'
        elif op == "--coreNum":
            coreNum = value
        elif op == "--hostFastaFileDir":
            hostFastaFileDir = value+'/'
    isExists = os.path.exists(outFileDir)
    if not isExists:
        os.mkdir(outFileDir)
        
    if virusFastaFileDir =='':
        print('please input the virusFastaFileDir.')
        return  
    elif not os.path.exists(virusFastaFileDir):
        print('please input a right virusFastaFileDir.')
        return        
    
    if hostFastaFileDir !='':
        if not os.path.exists(hostFastaFileDir):
            print('please input a right hostFastaFileDir.')
            return      
        hostKmerName = 'HostKmer_1to6_input'                    #inputHostFileDir
        hostKmerDir = outFileDir
        print('count kmer...')
        countKmer.getKmer(hostFastaFileDir,hostKmerDir,hostKmerName)                          #conutHostKmer
    else:
        print('please input the hostFastaFileDir.')
        return
        

                                                

    testVirusKmerDir = outFileDir                    #testVirusKmerDir
    
    
    countKmer.getKmer(virusFastaFileDir,outFileDir,testVirusKmerName)                        #conutVirusKmer
    dicVirusSeqLength = {}
    for dirVirus,b,virusNameList in os.walk(virusFastaFileDir):
        for eachVirusName in virusNameList:
            # print(eachVirusName)
            fileIn = open(dirVirus +'/'+ eachVirusName,'r')
            textIn = fileIn.readlines()
            fileIn.close()
            seq = ''
            for eachLine in textIn:
                if '>' not in eachLine:
                    seq += eachLine.strip('\n')
                    dicVirusSeqLength[eachVirusName] = len(seq)

    
    predictVirusHost(dicVirusSeqLength, testVirusKmerName, testVirusKmerDir, hostKmerDir, hostKmerName, outFileDir,int(coreNum))  #predictVirusHost
    
    

if __name__ == '__main__':
    main()
    #python3 phageHostPredictor.py --virusFastaFileDir='./exampleVirusGenome' --outFileDir='exampleOutput' --hostFastaFileDir='./exampleHostGenome'
    
    
    
    
    
    
    
    
    
    

