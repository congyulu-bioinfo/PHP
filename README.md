# Prokaryotic virus Host Predictor (PHP)
**Prokaryotic virus Host Predictor (PHP)** is a computational tool for host prediction of prokaryotic viruses based on **Gaussian Mixture Model (GMM)**. 
PHP takes the complete or partial genomic sequences of prokaryotic viruses as inputs. For each virus sequence, 
PHP automatically calculates the probability of host for 60,105 prokaryotic genomes, 
and takes the prokaryotic genome with the largest probability as the predicted host. 
PHP would output the name of prokaryotic genome with the largest score, and output the host score of all prokaryotic genomes. <br>  
Python packages `pandas`,`numpy`,`sklearn` are needed to be installed before Installation of phageHostPredictor. The program needs to run on the Linux operating system<br>  

Dependencies:
-----------
### Python 3, numpy, pandas, scikit-learn

    pip install pandas
    pip install numpy
    pip install -U scikit-learn


Usage:
-----------
Step 1: calculate the *K*-mer frequency of the host

    python3 countKmer.py --fastaFileDir  ./exampleHostGenome --kmerFileDir ./exampleOutput --kmerName HostKmer

Or use the simplify command

    python3 countKmer.py -f ./exampleHostGenome -d ./exampleOutput -n HostKmer

>`--fastaFileDir` or `-f`: The fasta file of prokaryotic genome sequences, one genome per file.<br>  
>`--kmerFileDir` or `-d`: The path of prokaryotic *K*-mer file.<br>  
>`--kmerName` or `-n`: The name of prokaryotic *K*-mer file.<br>  

Step 2: predict the infection relationship between the virus and the host

    python3 PHP.py --virusFastaFileDir ./exampleVirusGenome  --outFileDir ./exampleOutput  --bacteriaKmerDir ./exampleOutput  --bacteriaKmerName HostKmer

Or use the simplify command

    python3 PHP.py -v ./exampleVirusGenome  -o ./exampleOutput  -d ./exampleOutput  -n HostKmer

>`--virusFastaFileDir` or `-v`: The fasta file of query virus sequences, one virus genome per file.<br>  
>`--outFileDir` or `-o`: The path of temp files and result files.<br>  
>`--bacteriaKmerDir` or `-d`: The path of prokaryotic *K*-mer file.<br>  
>`--bacteriaKmerName` or `-n`: The name of prokaryotic *K*-mer file.<br>  


Interpretation of Result
-----------
After running the prediction program, you will see the output files Prediction_Maxhost and Prediction_Allhost in the output folder.<br>  
>In document `Prediction_Maxhost.csv`, The first column is the input virus，The third column is the highest score host and the second is the score of this host. <br>  
>In document `Prediction_Allhost.csv`, query viruses and scores for all prokaryotic genomes are given. <br>  


Users use their own data to customize the model
-----------
Users can use their own data to train customized models, `/exampleTrainingData/` provides sample data, each folder containing one pair of viruses and host genomes in exampleTrainingData. <br>
The virus needs to be saved in the `/phage/` folder, and the host needs to be saved in the `/host/` folder.<br>
Then run the following command to automatically train to get the model

    python3 PHP_UserTrain.py --trainDataDir ./exampleTrainingData --outFileDir ./exampleOutput --modelName UserModel.m

>`--trainDataDir` The training data folder, in which each pair of training data is saved in its own folder<br>  
>`--outFileDir` The path to save the trained model<br>  
>`--modelName` The name of the trained model<br>  

Users need to rename and replace the original PHP model to use the customized model<br>  <br>  <br>  




###predict short virus query contig
If your input sequences contain sequences shorter than 12500 bp, 
the program will automatically identify these short segments and use the corresponding model to predict.  
>if 12500bp< length，The sequence will be predicted using a full length model;<br>  
>if 7500bp < length <= 12500bp, The sequence will be predicted using a 10000bp model;<br>  
>if 4000bp < length <= 7500bp, The sequence will be predicted using a 5000bp model;<br>  
>if 2000 < length <= 4000bp, The sequence will be predicted using a 3000bp model;<br>  
>if length<=2000bp, The sequence will be predicted using a 1000bp model;<br>  






