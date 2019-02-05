#Necessary libraries

import numpy as np
import pandas as pd
import graphviz
import numexpr
import itertools
import argparse
import allFunctions

from subprocess import call
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.metrics import fbeta_score
from sklearn.metrics import accuracy_score
from allFunctions import randomForest
from allFunctions import rankInformative
from allFunctions import negativeOut
from allFunctions import binaryScore
from allFunctions import DT_cutOffs
from allFunctions import queryGenerator
from allFunctions import permutor
from allFunctions import fbetaTest


def runNSForest(rfTrees, threads, Median_Expression_Level, InformativeGenes, Genes_to_testing,betaValue):

    for column in dataDummy.columns[PrecolNum-1:PostcolNum]:

            ## Run Random Forest and get a ranked list
            Ranked_Features= randomForest(column, dataDummy, PrecolNum, rfTrees, threads)
            RankedList = rankInformative(Ranked_Features)

            ## Setup testArray for f-beta evaluation
            testArray = dataDummy[[column]]
            testArray.columns = ['y_true']

            #Rerank according to expression level and binary score
            Positive_RankedList_Complete = negativeOut(RankedList, column, medianValues, Median_Expression_Level)
            Binary_store_DF = pd.DataFrame()
            Binary_RankedList = binaryScore(Positive_RankedList_Complete, InformativeGenes, medianValues, column)

            Binary_score_store_DF_extra = Binary_store_DF.assign(clusterName = column)
            #print Binary_score_store_DF_extra
            Binary_score_store_DF = Binary_score_store_DF.append(Binary_score_store_DF_extra)


            #Get expression cutoffs for f-beta testing
            cut_dict = DT_cutOffs(Binary_RankedList,column)
            DT_cutoffs_store[column]=cut_dict

            #Generate expression queries and run those queries using fbetaTest() function
            queryInequalities=queryGenerator(Binary_RankedList, cut_dict)
            FullpermutationList = permutor(queryInequalities)
            #print len(FullpermutationList)
            f1_store = fbetaTest(FullpermutationList, column, testArray, betaValue)
            f1_store_1D.update(f1_store)



    #Report generation and cleanup

    f1_store_1D_df = pd.DataFrame() #F1 store gives all results.
    f1_store_1D_df = pd.DataFrame.from_dict(f1_store_1D, orient='index', dtype=str)
    f1_store_1D_df.columns = ["f-measure"]
    f1_store_1D_df['markerCount'] = f1_store_1D_df.index.str.count('&')
    f1_store_1D_df.reset_index(level=f1_store_1D_df.index.names, inplace=True)

    f1_store_1D_df_done= f1_store_1D_df['index'].apply(lambda x: pd.Series(x.split('&')))

    NSForest_Results_Table=f1_store_1D_df.join(f1_store_1D_df_done)

    NSForest_Results_Table_Fin = pd.DataFrame()
    NSForest_Results_Table_Fin = NSForest_Results_Table[NSForest_Results_Table.columns[0:4]]

    for i, col in enumerate(NSForest_Results_Table.columns[4:11]):
        splitResults= NSForest_Results_Table[col].astype(str).apply(lambda x: pd.Series(x.split('>=')))
        firstOnly = splitResults[0]
        Ascolumn = firstOnly.to_frame()
        Ascolumn.columns = [col]
        NSForest_Results_Table_Fin = NSForest_Results_Table_Fin.join(Ascolumn)


    NSForest_Results_Table_Fin.rename(columns={0:'clusterName'},inplace=True) #rename columns by position
    NSForest_Results_Table_Fin.sort_values(by=['clusterName','f-measure','markerCount'],ascending= [True, False, True], inplace = True)

    #Write outs
    Binary_score_store_DF.to_csv('Binary_scores_Supplmental_results.csv')
    NSForest_Results_Table_Fin.to_csv('NS-Forest_v2_results.csv')


    #Subsets of full results
    max_grouped = NSForest_Results_Table_Fin.groupby(by="clusterName")["f-measure"].max()
    max_grouped.df=pd.DataFrame(max_grouped)
    max_grouped.df.to_csv('NSForest_v2_maxF-scores.csv')

    NSForest_Results_Table_Fin["f-measureRank"] = NSForest_Results_Table_Fin.groupby(by="clusterName")["f-measure"].rank(ascending=False)
    topResults = NSForest_Results_Table_Fin["f-measureRank"] < 50
    NSForest_Results_Table_top = NSForest_Results_Table_Fin[topResults]
    NSForest_Results_Table_top.to_csv('NSForest_v2_topResults.csv')

#Pass parameters
parser = argparse.ArgumentParser(description='You can add a description here')
#parser.add_argument('-file',type=argparse.FileType('r'),help='filename')
parser.add_argument('-rfTrees',type=int,help='rfTrees number')
parser.add_argument('-threads',type=int,help='number of threads')
parser.add_argument('-Median_Expression_Level',type=int,help='Median_Expression_Level')
parser.add_argument('-InformativeGenes',type=int,help='InformativeGenes')
parser.add_argument('-Genes_to_testing',type=int,help='Genes_to_testing')
args = parser.parse_args()

print args.threads

dataFull = pd.read_table("Ab10k.tsv",index_col = 0)
#Creates dummy columns for one vs all Random Forest modeling
dataDummy = pd.get_dummies(dataFull, columns=["Clusters"], prefix = "", prefix_sep = "")

#Creates matrix of cluster median expression values
medianValues = dataFull.groupby(by="Clusters").median()
medianValues.to_csv('Function_medianValues.csv')

#Finding the number of clusters and printing that to screen (sanity check)
PrecolNum = len(dataFull.columns)
PostcolNum = len(dataDummy.columns)
adjustedColumns = PrecolNum-1
clusters2Loop=PostcolNum-PrecolNum

#Core analysis 
rankedDict =  {}  ###gives us the top ten features from RF
f1_store_1D = {}
Binary_score_store_DF=pd.DataFrame()
DT_cutoffs_store={}

####Random Forest parameters
#rfTrees=6 #Number of trees
#threads=1     #Number of threads to use, -1 is the greedy option where it will take all available CPUs/RAM

####Filtering and ranking of genes from random forest parameters

#InformativeGenes = 5 #How many top genes from the Random Forest ranked features will be evaluated for binariness
#Genes_to_testing = 6    #How many top genes ranked by binary score will be evaluated in permutations by fbeta-score (as the number increases the number of permutation rises exponentially!)

#### fbeta-score parameters
#Set values for fbeta weighting. 1 is default f-measure. close to zero is Precision, greater than 1 weights toward Recall

runNSForest(args.rfTrees, args.threads, args.Median_Expression_Level, args.InformativeGenes, args.Genes_to_testing,betaValue=0.5)
