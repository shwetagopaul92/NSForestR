#Necessary libraries

import numpy as np
import pandas as pd
import graphviz
import numexpr
import itertools
from subprocess import call

from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.metrics import fbeta_score
from sklearn.metrics import accuracy_score


#Function list

def randomForest(column, dataDummy,PrecolNum,rfTrees,threads):
    x_train = dataDummy[list(dataDummy.columns[0:PrecolNum-1])]
    names = dataDummy.columns[0:PrecolNum-1]
    y_train =dataDummy[column]
    rf = RandomForestClassifier(n_estimators = rfTrees, n_jobs = threads, random_state=123456)
    rf.fit(x_train, y_train)
    Ranked_Features=sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names),reverse=True)
    return Ranked_Features

def rankInformative(column,Ranked_Features):
    RankedList = []
    rankedDict = {}
    midcounter = 0
    for x in Ranked_Features:
        midcounter +=1
        RankedList.append(x[1])
        rankedDict[column] = RankedList
        if midcounter==30:
            break
    return RankedList

def negativeOut(x, column, medianValues,Median_Expression_Level):
    Positive_RankedList_Complete = []
    for i in x:
        if medianValues.loc[column, i] > Median_Expression_Level:
            print i
            print medianValues.loc[column, i]
            Positive_RankedList_Complete.append(i)
        else:
            print i
            print medianValues.loc[column, i]
            print "Is Right Out!"
    return Positive_RankedList_Complete

def binaryScore(Binary_store_DF,Genes_to_testing,Ranked_Features,Positive_RankedList_Complete, InformativeGenes, medianValues, column):
    Positive_RankedList=list(Positive_RankedList_Complete[0:InformativeGenes])
    Median_RF_Subset=medianValues.loc[:, Positive_RankedList]
    Rescaled_Matrix=pd.DataFrame()
    dataFull = pd.read_table("Ab10k.tsv",index_col = 0)
    dataDummy = pd.get_dummies(dataFull, columns=["Clusters"], prefix = "", prefix_sep = "")
    PrecolNum = len(dataFull.columns)
    PostcolNum = len(dataDummy.columns)
    clusters2Loop=PostcolNum-PrecolNum
    f1_store_1D = {}
    Binary_score_store_DF=pd.DataFrame()
    DT_cutoffs_store={}
    for i in Positive_RankedList:
        Target_value=medianValues.loc[column, i]
        Rescaled_values=Median_RF_Subset[[i]].divide(Target_value)
        Rescaled_Matrix=pd.concat([Rescaled_Matrix,Rescaled_values],axis=1)
    difference_matrix=Rescaled_Matrix.apply(lambda x: 1-x, axis=1)
    difference_matrix_clean = difference_matrix.where(difference_matrix > 0, 0)
    ColumnSums=difference_matrix_clean.sum(0)
    rescaled = ColumnSums/clusters2Loop

    # Double sort so that for ties, the RF ranking prevails!
    Ranked_Features_df=pd.DataFrame(Ranked_Features)
    Ranked_Features_df.rename(columns={1: 'Symbol'}, inplace=True)
    Ranked_Features_df_indexed=Ranked_Features_df.set_index("Symbol")
    rescaled_df=pd.DataFrame(rescaled)
    #global Binary_store_DF
    binaryAndinformation_Ranks=rescaled_df.join(Ranked_Features_df_indexed,lsuffix='_scaled', rsuffix='_informationGain')
    binaryAndinformation_Ranks.sort_values(by=['0_scaled','0_informationGain'],ascending= [False, False], inplace = True)
    Binary_ranked_Genes=binaryAndinformation_Ranks.index.tolist()
    Binary_RankedList=list(Binary_ranked_Genes[0:Genes_to_testing])
    Binary_scores=rescaled.to_dict()
    Binary_store_DF=Binary_store_DF.append(binaryAndinformation_Ranks)
    return Binary_RankedList

def DT_cutOffs(x, column):
    dataFull = pd.read_table("Ab10k.tsv",index_col = 0)
    dataDummy = pd.get_dummies(dataFull, columns=["Clusters"], prefix = "", prefix_sep = "")
    DT_cutoffs_store={}
    cut_dict = {}
    for i in x:
        filename=str(i)
        y_train =dataDummy[column]
        x_train = dataDummy[i]
        X = x_train[:, None]
        clf = tree.DecisionTreeClassifier(max_leaf_nodes=2)
        clf = clf.fit(X, y_train)
        threshold = clf.tree_.threshold
        cut_dict[i]=threshold[0]
    return cut_dict

def queryGenerator(x, cut_dict):
    queryList = []
    for i in x:
        str1 = i
        current_value = cut_dict.get(str1)
        queryString1 = str(str1)+'>='+ str(current_value)
        queryList.append(queryString1)
    return queryList

def permutor(x):
    binarylist2 = x
    combs = []
    for i in xrange(1, len(x)+1):
        els = [list(x) for x in itertools.permutations(binarylist2, i)]
        combs.extend(els)
    return combs

def fbetaTest(x, column,testArray, betaValue):
    fbeta_dict = {}
    dataFull = pd.read_table("Ab10k.tsv",index_col = 0)
    dataDummy = pd.get_dummies(dataFull, columns=["Clusters"], prefix = "", prefix_sep = "")
    for list in x:
        testArray['y_pred']= 0
        betaQuery = '&'.join(list)
        Ineq1=dataFull.query(betaQuery)
        testList=Ineq1.index.tolist()
        testArray.loc[testList, 'y_pred'] = 1
        f1 = fbeta_score(testArray['y_true'], testArray['y_pred'], average= 'binary', beta=betaValue)
        dictName = column+"&"+betaQuery
        fbeta_dict[dictName] = f1
    return fbeta_dict
