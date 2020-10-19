def userTrain(trainDataDir,outFileDir,modelName):
    import datetime,os,sys,getopt
    import pandas as pd
    import numpy as np
    from sklearn import mixture
    import joblib
    scriptPath = sys.path[0]+'/'

    
    
    for each in os.listdir(trainDataDir):
        kmerDir = trainDataDir+each+'/'
        virusKmerName = each+'.viruskmer'
        hostKmerName  = each+'.hostkmer'
        minusKmerName = each+'.minusKmer'
        #countVirusKmer
        virusDir = trainDataDir+each+'/phage/'
        if not os.path.exists(virusDir):
            print('There is no /phage/ folder in '+ trainDataDir+each+'/'+'\nbreak')
            return
        if len(os.listdir(virusDir)) == 0:
            print('There is no phage genome in '+ virusDir+'\nbreak')
            return
        
        os.system('python3 '+scriptPath+'countKmer.py --fastaFileDir  '+virusDir+' --kmerFileDir '+kmerDir+' --kmerName '+virusKmerName+' --coreNum -1')
        
        #countHostKmer
        hostDir = trainDataDir+each+'/host/'
        if not os.path.exists(hostDir):
            print('There is no /host/ folder in '+ trainDataDir+each+'/'+'\nbreak')
            return
        if len(os.listdir(hostDir)) == 0:
            print('There is no host genome in '+ hostDir+'\nbreak')
            return
        
        os.system('python3 '+scriptPath+'countKmer.py --fastaFileDir  '+hostDir+' --kmerFileDir '+kmerDir+' --kmerName '+hostKmerName+' --coreNum -1')


        #VirusKmerMinusHostKmer
        textVirus = pd.read_csv(kmerDir+virusKmerName, sep=",", header=None, index_col=0)
        textVirus.index = pd.Series([each])
        textHost  = pd.read_csv(kmerDir+hostKmerName, sep=",", header=None, index_col=0)
        textHost.index = pd.Series([each])
        minus = textVirus - textHost
        minus.to_csv(kmerDir+minusKmerName, index=True, header=False, sep=',', float_format='%.7f')
        
        
    #margeMinusKmerFile
    fileOut = open(outFileDir+'/trainingMinusKmer','w')
    for each in os.listdir(trainDataDir):
        kmerDir = trainDataDir+each+'/'
        for eachFile in os.listdir(kmerDir):
            if '.minusKmer' not in eachFile:
                continue
            fileIn = open(kmerDir+eachFile)
            textIn = fileIn.read()
            fileIn.close()
            fileOut.write(textIn)
    fileOut.close()
    print('kmer File Process completed.\nTraining Model...')
    
    #trainingModel
    text = pd.read_csv(outFileDir+'/trainingMinusKmer', sep=',', header=None, index_col=0).astype('float32')
    model = mixture.GaussianMixture(n_components=1,random_state=2)
    model.fit(text)
    joblib.dump(model, outFileDir + modelName)
    print('Model training completed.\n###  The model file '+modelName+' is saved in '+outFileDir)
    print('###  please move the model file to ./PHP/model/FullLength/ to use.')
    print('Users need to rename and replace the original PHP model to use the customized model\ndone')
    
    
def main():
    import pandas as pd
    import numpy  as np
    import datetime,os,sys,getopt
    trainDataDir = ''
    outFileDir = './exampleOutput/'
    modelName = 'UserModel.m'
    
    
    opts, args = getopt.getopt(sys.argv[1:], "ht:o:n:",["help=","trainDataDir=","outFileDir=","modelName="])
    for op, value in opts:
        if   op == "--trainDataDir"  or op == "-t":
            trainDataDir = value+'/'
        elif op == "--outFileDir"  or op == "-o":
            outFileDir = value+'/'
        elif op == "--modelName"  or op == "-n":
            modelName = value
        elif op == "--help"  or op == "-h":
            print('Users can use their own data to train customized models, /exampleTrainingData/ provides sample data, ')
            print('each folder containing one pair of viruses and host genomes in exampleTrainingData.')
            print('The virus needs to be saved in the /phage/ folder, and the host needs to be saved in the /host/ folder.')
            print('Then run the following command to automatically train to get the model')
            print()
            print('    python3 PHP_UserTrain.py --trainDataDir ./exampleTrainingData --outFileDir ./exampleOutput --modelName UserModel.m')
            print()
            print('--trainDataDir The training data folder, in which each pair of training data is saved in its own folder')
            print('--outFileDir The path to save the trained model')
            print('--modelName The name of the trained model')
            print()
            print('Users need to rename and replace the original PHP model to use the customized model')
            return
            
    if not os.path.exists(outFileDir):
        os.mkdir(outFileDir)
    if (not os.path.exists(trainDataDir)) or trainDataDir == '':
        print('please input a correct trainDataDir.\ndone.')
        return

    userTrain(trainDataDir,outFileDir,modelName)



if __name__ == '__main__':
    main()
    
    #python3 PHP_UserTrain.py --trainDataDir ./exampleTrainingData --outFileDir ./exampleOutput --modelName UserModel.m
    #python3 PHP_UserTrain.py -t ./exampleTrainingData -o ./exampleOutput -n UserModel.m







