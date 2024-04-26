#!/usr/bin/python
import decimal
import math
import os
import sys
import getopt
import pandas as pd
import multiprocessing
from numpy import log10
from scipy.stats import binomtest
from statsmodels.stats.proportion import proportion_confint
#import scipy #1.11.4
#print(scipy.__version__)

decimal.getcontext().prec = 400 # to change 50 to 300

''' this script programmed in python3 is designed for RMNE (random match Not exclude) to 
    calculate the probability to exclude one individual based on DNA alleles. The input file is from StrPhaser or VarSeqStitcher
    data format in the input file (frequency count): 
    CSF1PO	EAS	C,C,T,G,C,A,G,C,C,G;T,T,G;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT	17	T=A=S
    CSF1PO	EAS	C,C,T,A,C,G,G,C,C,C;T,T,GT;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT	12	T=A=S
        
    e.g.:  python AllelExclusionProb.py -d testData\stitched.MH.allele_count.Aalbers289.supPop.consisTag.txt -e testData\query.evidence.MH.allele.txt -s testData\query.suspect.MH.allele.txt -t 2 
    by Xuewen Wang
    Version:  Nov-28-2023
    Author: Xuewen Wang, Email: xwang.kib@gmail.com
'''
__author__="Xuewen Wang"

# to add clopper person for expected probability
#parse arguments
xapp=os.path.basename(__file__)
def usage():
    version="1.1.0"
    print(f"Usage: python {xapp} [options]")
    print(f"e.g., python3 {xapp} -d databaseAllele.txt -e evidenceAllele.txt -s suspectAllele.txt -t threads > outResults.txt")
    print("Options:")
    print("\t-h, --help: show this help message and exit")
    print("\t-d, --databaseAlle: string, required input file of dataBase with Allele count data, tab delimited")
    print("\t-e, --evidenceAlle: string, required input file of observed evidence Alleles, tab delimited")
    print("\t-s, --suspectAlle: string, required input file of suspect Alleles, tab delimited")
    print("\t-l, --log10: flag only, no value,output the probability after log10 conversion if -l is given.")
    print("\t-t, --threads: int, the number of parallelized computing threads, default 2")
    print("\t the result will be output to standard output/screen. use > to redirect to a file")
    print(f"Version: {version}, Dec,5th,2023")
    print("Support: xwang.kib@gmail.com")

def xpar(argv):
    #databaseAlle = "C:\\macrohaptype_testData\\testData\\databaseAllele.txt" #stitched.MH.allele_count.Aalbers289.supPop.consisTag.txt
    #evidenceAlle = "C:\\macrohaptype_testData\\testData\\evidenceAllele.txt" #query.evidence.MH.allele.txt
    #suspectAlle="C:\\macrohaptype_testData\\testData\\suspectAllele.txt" #query.suspect.MH.allele.txt
    #suspectAlle="C:\\macrohaptype_testData\\testData\\query.suspect.Excluded.MH.allele.txt"
    outfile = "C:\\macrohaptype_testData\\testData\\out_default.txt"
    num_threads = 2
    logarithm=False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:e:s:t:l", ["help","databaseAlle=", "evidenceAlle=", "suspectAlle=","threads=", "log10"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        exit(2)

    if not opts:
        # Print the usage message and exit
        print("Oops: No options provided")
        usage()
        sys.exit(0)


    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Tested to work well in  Python version: 3.9.13, scipy 1.11.4, statsmodels 1.16, numpy 1.21.0 and pandas version: 2.1.3")
            print("Your Python version:", sys.version)
            print("Your Pandas version:", pd.__version__)
            usage()
            MHdatabaseFormat='''
            input data format : 5 columns, tab separated
            CSF1PO	EAS	C,C,T,G,C,A,G,C,C,G;T,T,G;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT	17	T=A=S
            CSF1PO	EAS	C,C,T,A,C,G,G,C,C,C;T,T,GT;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT	12	T=A=S
            '''
            observedAlleFormat=    '''
            input data format :
            CSF1PO	C,C,T,G,C,A,G,C,C,G;T,T,G;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT
            CSF1PO	C,C,T,A,C,G,G,C,C,C;T,T,GT;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT
            '''
            suspectAlleFormat=    '''
            input data format : 3 columns of locus, allele, count(hom:2, het:1)
            CSF1PO	C,C,T,G,C,A,G,C,C,G;T,T,G;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT  1
            CSF1PO	C,C,T,A,C,G,G,C,C,C;T,T,GT;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT 1
            '''
            print("Database:\n"+MHdatabaseFormat, "\nObserved evidence:\n"+observedAlleFormat,"\nSuspected alleles:\n"+suspectAlleFormat, end="\n")
            exit(0)
        elif opt in ("-d", "--databaseAlle"):
            databaseAlle = arg
        elif opt in ("-e", "--sampleAlle"):
            evidenceAlle = arg
        elif opt in ("-s", "--suspectAlle"):
            suspectAlle = arg
        elif opt in ("-t", "--threads"):
            num_threads = arg
        elif opt in ("-l", "--log10"):
            logarithm = True
        else:
            assert False, "unexpected options"
    #print("#Input allele files: " +databaseAlle + " ; " +evidenceAlle+";"+suspectAlle+"\n")
    return databaseAlle, evidenceAlle, outfile, suspectAlle, num_threads, logarithm

def AlleleNormalziedFreq(df):
        ####calculate allele freq normalized to max 1 (but not equal to 1 or 0?)
        # input: dataframeOfMH, mkTotalAlleCountDict
        #CSF1PO  C,C,T,G,C,A,G,C,C,G,A,T,T,C,G,T,A,C,T,T,T,G,C,...  8
        sums = df.groupby(df.columns[0])[df.columns[2]].sum()
        # sums.reset_index() #to dataframe
        mkTotalAlleCountDict = sums.to_dict()
        mkPopAlleFreqDict = {}
        for x in df.index:
            mklocus = str(df.loc[x, 'Locus'])
            mkPopAlle = str(df.loc[x, 'Locus']) + ":"  + str(df.loc[x, 'Allele'])
            # CSF1PO:C,C,T,G,C,A,G,C,C
            count = df.loc[x, 'Count']
            mkPopAlleFreq = decimal.Decimal(count) / decimal.Decimal(mkTotalAlleCountDict.get(mklocus))
            mkPopAlleFreqDict[mkPopAlle] = mkPopAlleFreq
        return mkPopAlleFreqDict


def AlleleBionomialFreqHigh(df,confidLevel):
    # Clopper-Pearson binomial test to get the upbound of probability
    sums = df.groupby(df.columns[0])[df.columns[2]].sum()
    mkTotalAlleCountDict = sums.to_dict()
    mkPopAlleFreqDictBionom = {}
    nTotal=0
    kSee=0
    for x in df.index:
        mklocus = str(df.loc[x, 'Locus'])
        mkPopAlle = str(df.loc[x, 'Locus']) + ":" + str(df.loc[x, 'Allele'])
        # CSF1PO:C,C,T,G,C,A,G,C,C
        count = df.loc[x, 'Count']
        nTotal=int(mkTotalAlleCountDict.get(mklocus))
        kSee=int(count)
        # Calculate the probability of success
        #p = kSee/ nTotal

        confidInterval = binomtest(k=kSee, n=nTotal)  # hypothesized probability of 0.1
        ciHigh = confidInterval.proportion_ci(confidence_level=confidLevel).high
        #ciLow = confidInterval.proportion_ci(confidence_level=confidLevel).low
        ciPvalue=confidInterval.pvalue
        mkPopAlleFreqDictBionom[mkPopAlle] = ciHigh
        #reject if out of [cilow, cihigh]
    return mkPopAlleFreqDictBionom

def AlleleBionomialFreqHighStatsmodels(df,confidLevel):
    # Clopper-Pearson binomial test to get the upbound limit of probability via statsmodels
    sums = df.groupby(df.columns[0])[df.columns[2]].sum()
    mkTotalAlleCountDict = sums.to_dict()
    mkPopAlleFreqDictBionom = {}
    nTotal=0
    kSee=0
    for x in df.index:
        mklocus = str(df.loc[x, 'Locus'])
        mkPopAlle = str(df.loc[x, 'Locus']) + ":" + str(df.loc[x, 'Allele'])
        # CSF1PO:C,C,T,G,C,A,G,C,C
        count = df.loc[x, 'Count']
        nTotal=int(mkTotalAlleCountDict.get(mklocus))
        kSee=int(count)
        # Calculate the confidence interval
        confidInterval = proportion_confint(count=kSee, nobs=nTotal, alpha=1-confidLevel, method='normal') #'normal' is identical to scipy results; 'binom_test' is not working
       # print(confidInterval[1]) #(0.019585165831904984, 0.09981781924272187)
        ciHigh = confidInterval[1]
        #ciLow = confidInterval[0]
        mkPopAlleFreqDictBionom[mkPopAlle] = ciHigh
    return mkPopAlleFreqDictBionom

def suspectExcluded(sdf, obsAlleFreq):
    susUniqAlleles=[]
    for s in sdf.index:
         # CSF1PO	C,C,T,G,C,A,G,C,C,G,A,T,T,..;ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT 1
         susMkAlleSeq=sdf.iloc[s,0] + ":" + sdf.iloc[s,1]
         if obsAlleFreq.get(susMkAlleSeq) is None:
             susUniqAlleles.append(susMkAlleSeq.replace(":","\t"))
    if(len(susUniqAlleles) >=1):
         reportMessage = "The suspect is EXCLUDED because its unique allele(s) are not observed:\n "
         print(reportMessage)
         for sa in susUniqAlleles:
              print(sa)
         return True
    else:
         return False

def pExclusionAtOneLocus(alleList):
    # probability of exclusion (PE) at each locus: p=A1+A2+A3; q=1-p; PE=2pq+q2
    # 0.002, 0.1,0.4
    p=decimal.Decimal('0.0')
    for a in alleList:
        p += decimal.Decimal(a)
    q=decimal.Decimal('1')-p
    PE=2*p*q+q*q
    return PE

def pExclusionAtOneLocusLog(alleList):
    # probability of exclusion (PE) at each locus: p=A1+A2+A3; q=1-p; PE=2pq+q2
    # 0.002, 0.1,0.4
    p=decimal.Decimal('0.0')
    for a in alleList:
        p += decimal.Decimal(a)
    q=decimal.Decimal('1')-p
    #PE=2*p*q+q*q
    pq2=decimal(10) ** (log10(2)+log10(p)+log10(q))
    qq=decimal(10) **(2*log10(q))
    PE=pq2+qq
    return PE

def pExclusionAtOneLocusParalle(alleList):
    # probability of exclusion (PE) at each locus: p=A1+A2+A3; q=1-p; PE=2pq+q^2
    # [(CSF1PO	 [0.002, 0.1,0.4]),(vWA	 [0.2])]
    p=decimal.Decimal('0.0')
    mk=alleList[0]
    for a in alleList[1]:
        p += decimal.Decimal(a)
    if float(p) > 1.0:
        p = decimal.Decimal(1)
    q=decimal.Decimal('1')-p
    PE=2*p*q+q*q
    mkPE=(mk,PE)
    return mkPE # (CSF1PO, 0.99)

def pExclusionAtOneLocusParalleLog(alleList):
    # probability of exclusion (PE) at each locus: p=A1+A2+A3; q=1-p; PE=2pq+q2
    # [(CSF1PO	 [0.002, 0.1,0.4]),(vWA	 [0.2])]
    # added by AW: there is a bug; pq2 calculation gives Decimal runtime error sometimes.  we don't need logs at this stage, though I can't guarantee this is the right solution.
    return pExclusionAtOneLocusParalle(alleList)

    p=decimal.Decimal('0.0')
    mk=alleList[0]
    mkPE=(mk,0)
    for a in alleList[1]:
        p += decimal.Decimal(a)
        if float(p) > 1.0:
            p=decimal.Decimal(1)
    q=decimal.Decimal('1.0')-p
    #PE=2*p*q+q*q
    if float(p) > 0 and float(q) > 0:
        pq2=decimal.Decimal(10) **  (log10(decimal.Decimal(2))+log10(p)+log10(q))
        #pq2=decimal.Decimal(10) ** (decimal.Decimal(math.log10(2.0)) + decimal.Decimal(math.log10(float(p))) + decimal.Decimal(math.log10(float(q))))
        #qq=decimal.Decimal(10) **(decimal.Decimal(2)*log10(q))
        qq=decimal.Decimal(10) **(q.log10())
        PE=pq2+qq
        mkPE=(mk,PE)
    return mkPE # (CSF1PO, 0.99)

def pExclusionCombined(PEList):
    # combined probability of exclusion (CPE) for multiple loci n: CPE=1-[1-PE1]x(1-PE2)x...x(1-PEn)
    cpe=decimal.Decimal('1.0')
    cpelog=decimal.Decimal('1.0')
    for c in PEList:
        cpe *= (decimal.Decimal(1)-decimal.Decimal(c))
    cpe=decimal.Decimal(1)-cpe
    return cpe

def pExclusionCombinedAccuracy(PEList):
    # mathematical identity 1 - prod(Pi) = 1 - exp(sum(log(Pi)))
    cpe=decimal.Decimal('1.0')
    cpelog=decimal.Decimal('1.0')
    log_sum=decimal.Decimal(0)
    ctx = decimal.Context(prec=400)
    for c in PEList:
        # log_sum = sum(math.log(p) for p in P)
        ## cpe *= (decimal.Decimal(1)-decimal.Decimal(c)) equals the following
        pa=decimal.Decimal(1)-decimal.Decimal(c)
        # Calculate the sum of the logarithms
        #log_sum += math.log(pa) #error when p is 0
        if pa==0:
            continue
        log_sum += ctx.log10(pa) #from decimal package
    # Calculate 1 - exp(sum)
    #cpe = decimal.Decimal(1) - decimal.Decimal(math.exp(log_sum))
    cpe = decimal.Decimal(1) - ctx.exp(decimal.Decimal(log_sum))
    return cpe


def pExclusionCombinedLog(PEList):
    # combined probability of exclusion (CPE) for multiple loci n: CPE=1-[1-PE1]x(1-PE2)x...x(1-PEn)
    #use logarithms
    cpe=decimal.Decimal('1.0')
    cpelog=decimal.Decimal('1.0')
    product=decimal.Decimal('1.0')
    for c in PEList:
        diff=decimal.Decimal(1.0)-decimal.Decimal(c)
        if diff==0:
            continue
        product += log10(diff)
        #convert back
        cpe=decimal.Decimal(10.0)**product
        cpe=decimal.Decimal(1.0)-cpe
    return cpe

def getMHdataBaseSummary(df):
    print("#### Database summary: ")
    print("#Data format looks like:")
    print(df.head())
    popGroups=df[df.columns[1]].nunique()
    print("#Total populations:\t", popGroups)
    locusGroups=df[df.columns[0]].nunique()
    print("#Total MH loci:\t", locusGroups)
    TotalMH = df[df.columns[3]].sum()
    print("#Total MH alleles:\t", TotalMH, "\n")

def getMHallesSummary(df):
    print("#Data format looks like:")
    print(df.head())
    locusGroups=df[0].nunique()
    print("#Total  loci:\t", locusGroups)
    TotalMH = df.shape[0]
    print("#Total alleles:\t", TotalMH, "\n")



if __name__ == "__main__":
    popAlleCount, evidenceAlle, outfile,suspectAlle,num_threads,logarithm = xpar(sys.argv[1:])
    mkAlleCountDict = {}
    mkAlleFreqinOnePop={}
    mkAlleFreqinOnePopBinom={}
    obsAlleProbabilityBionomStat={}
    obsAlleFreq={}
    susAlleFreq={}
    num_threads=int(num_threads)
    confidLevel=0.95


    #prepare outfile
    '''
    if outfile != "out_default.txt":
        OUTH = open(outfile, "w")
        #OUTH.write("#")
    '''

    #### read MH database file to dataframe
    df = pd.read_table(popAlleCount, header=None)
    df.dropna(inplace=True)
    df.columns = ['Locus', 'Population', 'Allele', 'Count', 'Consistency']
    '''
          Locus   Population   Allele      Count Consistency
    0     CSF1PO  EAS          C,C,T,G...  17    T=A=S
    1     CSF1PO  EAS          C,C,T,A...  12    T=A=S
    '''
    # summary of MH database
    #getMHdataBaseSummary(df)
    popGroups = df.groupby(df.columns[1])  # by population, dict: pop:rowIndex1,rowIndex2

    #### read observed evidence file to dataframe
    #print("#### Observed evidence data: ")
    odf=pd.read_table(evidenceAlle, header=None)
    #getMHallesSummary(odf)
    obsAlleFreq={}
    for o in odf.index:
        obsAlleFreq[odf.iloc[o,0] + ":" + odf.iloc[o,1]] = 1 # 'CSF1PO:C,C,T,A,C,G,G,C' : 1

    #### read suspect allelic file
    sdf = pd.read_table(suspectAlle, header=None)
    #print("#### Suspect data: ")
    #getMHallesSummary(sdf)
    sdf.columns=['Locus','Allele','Count']
    #sdf['Count']=1
    for s in sdf.index:
        susAlleFreq[sdf.iloc[s,1] + ":" + sdf.iloc[s,1]] = sdf.iloc[s,2] # 'CSF1PO:C,C,T,A,C,G,G,C' : 1

    #### suspect alleles are not present or the probability in observed alleles
    #print("#### Report")
    if suspectExcluded(sdf, obsAlleFreq):
        #report and exit
        sys.exit(0)
    else:
        #print("# Probability of Exclusion (PE) / Inclusion (PI) of the suspect:")
        print("Population\tLocus" + "\tProbability_of_Exclusion(PE)" + "\tProbability_of_Inclusion(PI)"+"\tAdjusted_PE\u00b2\tAdjusted_PI\u00b2"+"\tAdjusted_PE\u00b3\tAdjusted_PI\u00b3")

        ## get exlcusion prob in each population in MH database
        for pop, pgroup in popGroups:
            #print("Population:\t" + pop)
            '''
                 Locus Population  Allele                                     Count  Consistency
            0    CSF1PO  SAS  C,C,T,G,C,A,G,C,C,G,A,T,T,C,G,T,A,C,T,T,T,G,C,...  6  T=A=S
            1    CSF1PO  SAS  C,A,A,G,C,A,A,T,C,C,A,C,G,C,A,T,C,C,T,T,C,G,C,...  4  T=A=S
            '''

            popdf = pgroup.reset_index(drop=True) # population dataframe with new index
            popdf=popdf.drop(popdf.columns[[1,4]], axis=1)
            '''    Locus                                             Allele  Count
                0  CSF1PO  C,C,T,G,C,A,G,C,C,G,A,T,T,C,G,T,A,C,T,T,T,G,C,...      6
                1  CSF1PO  C,A,A,G,C,A,A,T,C,C,A,C,G,C,A,T,C,C,T,T,C,G,C,...      4
            '''

            # update MH database: adding new alleles or update the count
            popdfUpdated = pd.merge(popdf, sdf, on=['Locus', 'Allele'], how='outer', suffixes=('_df1', '_df2'))
            popdfUpdated['Count'] = popdfUpdated['Count_df1'].fillna(0) + popdfUpdated['Count_df2'].fillna(0)
            popdfUpdated.drop(['Count_df1', 'Count_df2'], axis=1, inplace=True)
            popdfUpdated.reset_index()
            # Locus   Allele                                             Count
            # CSF1PO  C,C,T,G,C,A,G,C,C,G,A,T,T,C,G,T,A,C,T,T,T,G,C,...  17
            '''for popu in popdfUpdated.index:
                #checkalle="G,G,T,G,A,A,A,T,C,A,C,T,T,T,C,A,T,A,T,G,G,A,C,C,C,A,C,T,T,C,C,T,T,G,A,C,C,T,T,A,A,A,G,T,G,A,A,G,T,C,G,G,C,G,C,C,T,G,C,G,T,G,C,C,G,G,C,T,C,A,C,G,T,G,C,G,T,G;T,G,AAC,CA,T;TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATTAGATAGACAGACAGACAGATAGA"
                checkalle="G,G,T,G,A,T,T,T,C,A,C,T,T,T,C,A,T,A,T,G,G,A,C,C,C,A,C,T,T,C,C,T,T,G,A,C,C,T,T,A,A,A,G,T,G,A,A,G,T,C,G,G,C,G,C,C,T,G,C,G,T,G,C,C,G,G,C,T,C,A,C,G,T,G,C,G,T,G;T,G,AAC,CA,T;TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATTAGATAGACAGACAGACAGATAGA"
                if popdfUpdated.iloc[popu, 0]+":"+popdfUpdated.iloc[popu, 1]== 'vWA'+":"+ checkalle :
                    print(popdfUpdated.loc[popu].to_list)
            '''
            #calculate allele freq ratio normalized to max 1
            mkAlleFreqinOnePop = AlleleNormalziedFreq(popdfUpdated)
            mkAlleFreqinOnePopBinom=AlleleBionomialFreqHigh(popdfUpdated, confidLevel)
            mkAlleFreqinOnePopBinomStat=AlleleBionomialFreqHighStatsmodels(popdfUpdated, confidLevel)

            # get the observed evidence AlleSeq frequency in Database or not, if not: p=0
            obsAlleProbability={}
            obsAlleProbabilityBionom={}
            popMkAlleSeq_probBionom=0.0
            popMkAlleSeq_probBionomStat=0.0
            for obsAlle in obsAlleFreq:
                mkname=obsAlle.split(':')[0]
                if mkAlleFreqinOnePop.get(obsAlle) is None:
                    popMkAlleSeq_prob = 0
                else:
                    popMkAlleSeq_prob = mkAlleFreqinOnePop.get(obsAlle) #0.01492537313432835820895522388
                    popMkAlleSeq_probBionom=mkAlleFreqinOnePopBinom.get(obsAlle)
                    popMkAlleSeq_probBionomStat=mkAlleFreqinOnePopBinom.get(obsAlle)

                if mkname in obsAlleProbability:
                        obsAlleProbability[mkname].append(popMkAlleSeq_prob)
                        obsAlleProbabilityBionom[mkname].append(popMkAlleSeq_probBionom)
                        obsAlleProbabilityBionomStat[mkname].append(popMkAlleSeq_probBionomStat)
                else:
                    probList=[popMkAlleSeq_prob]
                    obsAlleProbability[mkname]=probList
                    probListBionom = [popMkAlleSeq_probBionom]
                    obsAlleProbabilityBionom[mkname] = probListBionom
                    probListBionomStat = [popMkAlleSeq_probBionomStat]
                    obsAlleProbabilityBionomStat[mkname] = probListBionomStat

            # prepare data to calculate prob at a locus for all evidence loci
            PEList=[]
            mk_alleProbList=[]
            PEListBionom = []
            PEListBionomStat = []
            mk_alleProbListBionom = []
            mk_alleProbListBionomStat = []

            for mkname in obsAlleProbability:
                alleProbList=obsAlleProbability.get(mkname)
                alleProbListBionom = obsAlleProbabilityBionom.get(mkname)
                alleProbListBionomStat = obsAlleProbabilityBionomStat.get(mkname)

                ## single thread
                #PE =pExclusionAtOneLocus(alleProbList) #probability of exclusion (PE)
                #print(mkname+"\t",PE)
                #PEList.append(PE)
                mk_alleProbList.append((mkname,alleProbList)) # (CSF1PO,[0.002,0.1,0.4]),(vWA, [0.2])
                mk_alleProbListBionom.append((mkname, alleProbListBionom))
                mk_alleProbListBionomStat.append((mkname, alleProbListBionomStat))

            #paralleling calculate PEs for markers in the same population
            pool = multiprocessing.Pool(num_threads)
            #locusPEList=pool.map(pExclusionAtOneLocusParalle, mk_alleProbList)
            locusPEList = pool.map(pExclusionAtOneLocusParalleLog, mk_alleProbList)
            locusPEListBionom = pool.map(pExclusionAtOneLocusParalleLog, mk_alleProbListBionom)
            locusPEListBionomStat = pool.map(pExclusionAtOneLocusParalleLog, mk_alleProbListBionomStat)
            # [('CSF1PO', Decimal('0.9986282578875171467764060357')), ('vWA', Decimal('0.9999126561271726788365796140'))]
            pool.close()
            pool.join()

            #### output
            ## locus Probe list: outProbList
            for m in range(0,len(locusPEList)):
                mk=locusPEList[m][0]
                PE=locusPEList[m][1]
                PEBionom = locusPEListBionom[m][1]
                PEBionomStat = locusPEListBionomStat[m][1]
                PE1minus=decimal.Decimal(1) - PE
                PE1minusBionom = decimal.Decimal(1) - PEBionom
                PE1minusBionomStat = decimal.Decimal(1) - PEBionomStat

                if logarithm:
                    #log10 output
                    PEplus = pop + "\t" + mk + "\t" + str(log10(PE)) + "\t" + str(log10(PE1minus))+"\t"+str(log10(PEBionom))+"\t" + str(log10(PE1minusBionom))+"\t"+str(log10(PEBionomStat))+"\t" + str(log10(PE1minusBionomStat))
                    print(PEplus)
                else:
                    #normal output: default
                    #PEplus = pop + "\t" + mk + "\t" + str(PE) + "\t" + str(PE1minus)+"\t"+str(PEBionom)+"\t" + str(PE1minusBionom)+"\t"+str(PEBionomStat)+"\t" + str(PE1minusBionomStat)
                    print(f"{pop}\t{mk}\t{PE:.20E}\t{PE1minus:.20E}\t{PEBionom:.20E}\t{PE1minusBionom:.20E}\t{PEBionomStat:.20E}\t{PE1minusBionomStat:.20E}")

                #print(PEplus)
                PEList.append(PE)
                PEListBionom.append(PEBionom)
                PEListBionomStat.append(PEBionomStat)

            ##combined probability of exclusion (CPE)
            CPE=0.0
            #CPE=pExclusionCombined(PEList)
            CPE=pExclusionCombinedAccuracy(PEList) #updated
            CPEBionom=pExclusionCombined(PEListBionom)
            CPEBionomStat = pExclusionCombined(PEListBionomStat)
            #CPE = pExclusionCombinedLog(PEList)
            #CPEBionom=pExclusionCombinedLog(PEListBionom)
            #CPEBionomStat = pExclusionCombinedLog(PEListBionomStat)

            CPI=decimal.Decimal(1)-decimal.Decimal(CPE) #Decimal('1')
            CPIBionom = decimal.Decimal(1) - decimal.Decimal(CPEBionom)
            CPIBionomStat = decimal.Decimal(1) - decimal.Decimal(CPEBionomStat)
            if logarithm:
                CpPlus=print(f"{pop}\tAll_loci\t{log10(CPE):.20E}\t{log10(CPI):.20E}\t{log10(CPEBionom):.20E}\t{log10(CPIBionom):.20E}\t{log10(CPEBionomStat):.20E}\t{log10(CPIBionomStat):.20E}")
            else:
                print(f"{pop}\tAll_loci\t{CPE:.20E}\t{CPI:.20E}\t{CPEBionom:.20E}\t{CPIBionom:.20E}\t{CPEBionomStat:.20E}\t{CPIBionomStat:.20E}") #{CPI:.28E}



   #OUTH.close()
