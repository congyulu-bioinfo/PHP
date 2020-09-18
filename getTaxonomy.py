def getTop100Host():
    fileIn = open('./virusKmer_Prediction_Allhost','r')
    textIn = fileIn.readlines()
    fileIn.close()
    hostNameList = textIn[0].strip('\n,').split(',')[1:]
    print('bactNum:',len(hostNameList))
    fileOut = open('./List_consensus','w')
    for eachLine in textIn[1:]:
        scoreList = eachLine.strip('\n,').split(',')[1:]
        virusName = eachLine.strip('\n,').split(',')[0]
        print(virusName)
        
        LIST = list(zip(hostNameList,scoreList))
        sortedList = sorted(LIST,key=lambda x:float(x[1]),reverse=True)
        fileOut.write(virusName)
        for each in sortedList[:100]:
            fileOut.write('\t'+str([each[0],each[1]]))
        fileOut.write('\n')
    fileOut.close()

def TransID2GenusID():
    fileIn = open('./60105ID_allRank_Complete_addGCA','r')
    textIn = fileIn.readlines()
    fileIn.close()
    dicID2Genus = {}
    dicGenus2ID = {}
    for eachLine in textIn:
        eachLine = eachLine.strip('\n')
        norankID = eachLine.split('\t')[0].split(',')[0].split('.')[0]
        
        if RANK == 'Genus':
            rankID = eachLine.split('\t')[-2]
        if RANK == 'Family':
            rankID = eachLine.split('\t')[-3]
        if RANK == 'Order':
            rankID = eachLine.split('\t')[-4]
        if RANK == 'Class':
            rankID = eachLine.split('\t')[-5]
        if RANK == 'Phylum':
            rankID = eachLine.split('\t')[-6]
        if RANK == 'Domain':
            rankID = eachLine.split('\t')[-7]
        dicID2Genus[norankID] = rankID
        dicGenus2ID[rankID] = norankID
    
    fileIn = open('./List_consensus','r')
    textIn = fileIn.readlines()
    fileIn.close()
    fileOut = open('./List_consensus_Genus','w')
    for eachLine in textIn:
        virusName = eachLine.split('\t')[0]
        print(virusName)
        fileOut.write(virusName)
        for each in eachLine.split('\t')[1:]:
            print(each)
            hostNorankID = eval(each)[0]
            score = eval(each)[1]
            IDTransed = dicGenus2ID[dicID2Genus[hostNorankID]]
            fileOut.write('\t'+str([IDTransed,str(score)]))
        fileOut.write('\n')
    fileOut.close()

def makeConsensus(NUM):
    from collections import Counter
    fileIn = open('./List_consensus_Genus','r')
    textIn = fileIn.readlines()
    fileIn.close()
    
    fileOut = open('./List_consensus_Genus_top','w')
    for eachLine in textIn:
        virusName = eachLine.split('\t')[0]
        LIST = []
        for each in eachLine.split('\t')[1:][:NUM]:
            LIST.append(eval(each)[0])
        
        collection_words = Counter(LIST)
        most_counterNum = collection_words.most_common(1)
        fileOut.write(virusName+'\t'+str(most_counterNum[0][1])+'\t'+str(most_counterNum[0][0])+'\n')
    fileOut.close()



if __name__ == '__main__':
    import pandas as pd
    import numpy  as np
    getTop100Host()
    for RANK in ['Genus','Family','Order','Class','Phylum','Domain']:
        TransID2GenusID()
        NUM = 5
        makeConsensus(NUM=NUM)

















