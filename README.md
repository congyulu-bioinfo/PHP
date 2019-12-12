# phageHostPredictor
PhageHostPredictor (PHP) is a computational tool for host prediction of phages. PHP takes the complete or partial genomic sequences of phages as inputs. For each phage, PHP automatically calculates the probability of host for 36,600 bacterial species, and takes the bacterial species with the largest probability as the predicted host. PHP would output the name, the probability, and output the host probability of all bacterial genomes. <br>  
Python packages "pandas","numpy","sklearn" are needed to be installed before Installation of phageHostPredictor. The program needs to run on the Linux operating system<br>  


Dependencies:
-----------
### Python 3, numpy, pandas, scikit-learn

    pip install pandas
    pip install numpy
    pip install -U scikit-learn


Usage:
-----------
### decompression the model file
    
    for tar in ./model/*.tar.gz; do tar xvf $tar -C ./model; done
    
### example command
If you want to predict the infection relationship between the virus you provide and the host

    python3 phageHostPredictor.py --virusFastaFileDir='./exampleVirusGenome' --outFileDir='./exampleOutput'  --hostFastaFileDir='./exampleHostGenome'
    
--virusFastaFileDir	The fasta file of query virus sequences, one virus genome per file.<br>  
--outFileDir	The dir of temp files and result files.<br>  
--hostFastaFileDir	The fasta file of host sequences, one host genome per file.<br>  


### simplified example command, 
you can omit the output path. The result file will be stored at ./exampleOutput and and one process will be used.

    
    python3 phageHostPredictor.py --virusFastaFileDir='./exampleVirusGenome' --hostFastaFileDir='./exampleHostGenome'


    
    
Interpretation of Result
-----------
After running the prediction program, you will see the output files Prediction_Maxhost and Prediction_Allhost in the output folder.

In document Prediction_Maxhost, The first column is the input virus，The third column is the highest probability host and the second is the probability of this host. 

In document Prediction_Allhost, query viruses and probability for all bacterial genomes are given. 


Short virus query sequence
-----------
If your input sequences contain sequences shorter than 12500 bp, the program will automatically identify these short segments and use the corresponding model to predict.  
if 12500<length，The sequence will be predicted using a full length model;<br>  
if 7500<length<=12500, The sequence will be predicted using a 10000bp model;<br>  
if 4000<length<=7500, The sequence will be predicted using a 5000bp model;<br>  
if 2500<length<=4000, The sequence will be predicted using a 3000bp model;<br>  
if 1500<length<=2500, The sequence will be predicted using a 2000bp model;<br>  
if 750<length<=1500, The sequence will be predicted using a 2000bp model;<br>  
if length<=750, The sequence will be predicted using a 500bp model;<br>  

