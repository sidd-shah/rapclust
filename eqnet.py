from __future__ import print_function
import time
def buildNetFile(sampdirs, netfile, cutoff, writecomponents=False):
    import itertools
    import pandas as pd
    import numpy as np
    import os

    sep = os.path.sep
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]

    quant = None
    for sffile in sffiles:
        if quant is None:
            quant = pd.read_table(sffile)
            quant.set_index('Name', inplace=True)
        else:
            quant2 = pd.read_table(sffile)
            quant2.set_index('Name', inplace=True)
            quant += quant2

    tnames = []
    weightDict = {}
    diagCounts = np.zeros(len(quant['TPM'].values))

    tot = 0
    eqfiles = [sep.join([sd, 'aux/eq_classes.txt']) for sd in sampdirs]

    firstSamp = True
    numSamp = 0
    eqClasses = {}
    for eqfile in eqfiles:
        with open(eqfile) as ifile:
            numSamp += 1
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
            if firstSamp:
                for i in xrange(numTran):
                    tnames.append(ifile.readline().rstrip())
            else:
                for i in xrange(numTran):
                    ifile.readline()

            for i in xrange(numEq):
                toks = map(int, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:-1])
                count = toks[-1]
                if tids in eqClasses:
                    eqClasses[tids] += count
                else:
                    eqClasses[tids] = count

            firstSamp = False

    tpm = quant.loc[tnames, 'TPM'].values / numSamp
    estCount = quant.loc[tnames, 'NumReads'].values
    efflens = quant.loc[tnames, 'EffectiveLength'].values
    epsilon =  np.finfo(float).eps
    for tids, count in eqClasses.iteritems():
        denom = sum([tpm[t] for t in tids])
        tot += count
        for t1, t2 in itertools.combinations(tids,2):
            #tpm1 = tpm[t1]
            #tpm2 = tpm[t2]
            #w = count * ((tpm1 + tpm2) / denom)
            if (t1, t2) in weightDict:
                weightDict[(t1, t2)] += count
            else:
                weightDict[(t1, t2)] = count
        for t in tids:
            #if (estCount[t] <= cutoff):
            #    continue
            #diagCounts[t] += count * (tpm[t] / denom)
            diagCounts[t] += count


    print("total reads = {}".format(tot))
    maxWeight = 0.0
    prior = 0.1
    edgesToRemove = []
    for k,v in weightDict.iteritems():
        c0, c1 = diagCounts[k[0]], diagCounts[k[1]]
        #w = (v + prior) / (min(c0, c1) + prior)
        if c0 + c1 > epsilon and c0 > cutoff and c1 > cutoff:
            w = v / min(c0, c1)
            weightDict[k] = w
            if w > maxWeight:
                maxWeight = w
        else:
            edgesToRemove.append(k)

    for e in edgesToRemove:
        del weightDict[e]

    tnamesFilt = []
    relabel = {}
    for i in xrange(len(estCount)):
        if (diagCounts[i] > cutoff):
            relabel[i] = len(tnamesFilt)
            tnamesFilt.append(tnames[i])
            weightDict[(i, i)] = 1.1

    import networkx as nx
    G = nx.Graph() if writecomponents else None
    with open(netfile, 'w') as ofile:
        writeEdgeList(weightDict, tnames, ofile, G)

    if G is not None:
        clustFile = netfile.split('.net')[0] + '.clust'
        print("Writing connected components as clusters to {}".format(clustFile))
        with open(clustFile, 'w') as ofile:
            cc = nx.connected_component_subgraphs(G)
            for c in cc:
                ofile.write('{}\n'.format('\t'.join(c.nodes())))

def writeEdgeList(weightDict, tnames, ofile, G):
    useGraph = G is not None
    for k,v in weightDict.iteritems():
        ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
        if useGraph:
            G.add_edge(tnames[k[0]], tnames[k[1]])


def writePajek(weightDict, tnames, relabel, ofile):
    with open(netfile, 'w') as ofile:
        ofile.write("*Vertices\t{}\n".format(len(tnamesFilt)))
        for i, n in enumerate(tnamesFilt):
            ofile.write("{}\t\"{}\"\n".format(i, n))
        ofile.write("*Edges\n")
        print("There are {} edges\n".format(len(weightDict)))
        for k,v in weightDict.iteritems():
            ofile.write("{}\t{}\t{}\n".format(relabel[k[0]], relabel[k[1]], v))
            #ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
            #if k[0] != k[1]:
            #    ofile.write("{}\t{}\t{}\n".format(tnames[k[1]], tnames[k[0]], v))

class EquivCollection(object):
    def __init__(self):
        self.tnames = []
        self.eqClasses = {}
        self.hasNames = False

    def setNames(self, names):
        self.tnames = names
        self.hasNames = True

    def add(self, tids, count):
        if tids in self.eqClasses:
            self.eqClasses[tids] += count
        else:
            self.eqClasses[tids] = count

def readEqClass(eqfile, eqCollection):
    with open(eqfile) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
        if not eqCollection.hasNames:
            tnames = []
            for i in xrange(numTran):
                tnames.append(ifile.readline().rstrip())
            eqCollection.setNames(tnames)
        else:
            for i in xrange(numTran):
                ifile.readline()

        for i in xrange(numEq):
            toks = map(int, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = tuple(toks[1:-1])
            count = toks[-1]
            eqCollection.add(tids, count)

def getCountsFromEquiv(eqCollection):
    countDict = {}
    tn = eqCollection.tnames
    for tids, count in eqCollection.eqClasses.iteritems():
        for t in tids:
            if tn[t] in countDict:
                countDict[tn[t]] += count
            else:
                countDict[tn[t]] = count
    # ensure no division by 0
    for t in eqCollection.tnames:
        if t in countDict:
            countDict[t] += 1.0
        else:
            countDict[t] = 1.0
    return countDict

def flattenClusters(infile, outfile):
    with open(outfile, 'w') as ofile:
        with open(infile) as ifile:
            for i,l in enumerate(ifile):
                toks = l.rstrip().split()
                cname = "cluster{}".format(i)
                for t in toks:
                    ofile.write("{}\t{}\n".format(cname, t))

def calculate_prob(edge, ambigCounts, eqCollection):
    num = edge['count']
    count1 = ambigCounts[edge['t1']]
    count2 = ambigCounts[edge['t2']]
    denom = count1 + count2 - num
    edge['prob'] = float(num)/float(denom)
    edge['total'] = denom
    edge['pre'] = 2* edge['prob']*edge['prob'] * (1- edge['prob'])
    return edge

def buildEdgeProbTable(eqCollection):
    binom_classes = {}
    count = 1
    for row in eqCollection.eqClasses.iteritems():
        if count < 0:
            break
        else:
            tup = row[0]
            read_count = row[1]
            if(len(tup)>=2):
                for i1,contig1 in enumerate(tup):
                    for i2,contig2 in enumerate(tup):
                        if(contig1!=contig2) and i2>i1:
                            current_count = read_count
                            if (contig1, contig2) in binom_classes:
                                current_count = binom_classes[(contig1, contig2)] + read_count
                                binom_classes[(contig1, contig2)] = current_count
                            elif (contig2, contig1) in binom_classes:
                                current_count = binom_classes[(contig2, contig1)] + read_count
                                binom_classes[(contig2, contig1)] = current_count
                            else:
                                 binom_classes[(contig1, contig2)] = current_count
                                
            count +=1
    return binom_classes

def generate_dict(master):
    mast_dict = {}
    for i,row in master.iterrows():

        row_dict ={}
        key = (row['t1'], row['t2'])
        if key in mast_dict:
            row_dict = mast_dict[key]        
            row_row_dict = {}
            if row['Condition'] in row_dict:
                row_row_dict = row_dict[row['Condition']]
            row_row_dict[row['Folder']] = {
                'count' : row['count'],
                'total': row['total'],
                'prob':row['prob'],
                'pre':row['pre']
            }
            row_dict[row['Condition']] = row_row_dict
        else:
            row_dict[row['Condition']] = {}
            row_dict[row['Condition']][row['Folder']] = {}
            folder = row_dict[row['Condition']][row['Folder']]
            folder['count'] = row['count']
            folder['total']= row['total']
            folder['prob']=row['prob']
            folder['pre'] =row['pre']
            mast_dict[key] = row_dict
    return mast_dict

def write_to_file(list, file):
    for line1 in list:
        for line2 in line1:
            file.write(','.join(map(str, line2)))
            file.write("\n")

def filterGraph(expDict, netfile, ofile):
    import os
    import pandas as pd
    import math

    # Get just the set of condition names
    conditions = expDict.keys()
    print("conditions = {}".format(conditions))

    #for cond in conditions:
    #    sailfish[cond] = collections.defaultdict(float)
    #    for sample in samples:
    #        fileDir = cond+sample
    #        filePath = os.path.sep.join([samplesPath, fileDir, "quant.sf"])
    #        print(filePath)
    #        with open(filePath) as f:
    #            data = pd.read_table(f, header=0)
    #            for i in range(len(data)):
    #                name = data['Name'][i]
    #                numreads = data['NumReads'][i]
    #                sailfish[cond][name] += numreads
    #    # To avoid divide by 0 error
    #    for name in sailfish[cond]:
    #        sailfish[cond][name] += 1
    eqClasses = {}
    masterDf = None
    try:
        masterDf = pd.read_csv('/home/ec2-user/proj1/sailfish_quant/masterdf.csv')
    except:
        pass
    if masterDf is None:
        for cond in conditions:
            print("Expression Dict",expDict[cond])
            for sampNum, sampPath in expDict[cond].iteritems():
                print("Samp Num {0}, Samp Path {1}".format(sampNum, sampPath))
                if cond not in eqClasses:
                    eqClasses[cond] = EquivCollection()
                eqPath = os.path.sep.join([sampPath, "aux", "eq_classes.txt"])
                print("Eq Path: ", eqPath)
                readEqClass(eqPath, eqClasses[cond])
                print("\nPRINTING EQ CLASSES ", len(eqClasses[cond].eqClasses))
                myEqCollection = EquivCollection()
                readEqClass(eqPath, myEqCollection)
                counts_dict = buildEdgeProbTable(myEqCollection)
                local_ambig = getCountsFromEquiv(myEqCollection)
                df = pd.DataFrame(columns=['Condition','Folder','t1','t2','count','total', 'prob','pre'])
                df['t1'] = [myEqCollection.tnames[item[0]] for item in counts_dict.keys()]
                df['t2'] = [myEqCollection.tnames[item[1]] for item in counts_dict.keys()]
                df['count'] = counts_dict.values()
                df['Condition'] = cond
                df['Folder'] = sampPath
                df = df.apply(calculate_prob,  axis=1,args=([local_ambig, myEqCollection]))
                masterDf = pd.concat([masterDf, df])
        ambigCounts = {cond : getCountsFromEquiv(eqClasses[cond]) for cond in conditions}

        sailfish = {}
        for cond in conditions:
            sailfish[cond] = ambigCounts[cond]         
        masterDf.to_csv('masterdf.csv')
    mast_dict = generate_dict(masterDf)
    print("Master DF Length", len(masterDf))
    print ("Done Reading")
    count = 0
    numTrimmed = 0
    print("NETFILE ", netfile)
    print("OFILE ", ofile)
    import csv   
    big_file = open('biggest_file.csv', 'wb')
    big_df = pd.DataFrame()
    big_list = []
    wr =csv.writer(big_file)
    with open(netfile) as f, open(ofile, 'w') as ofile:
        data  = pd.read_table(f, header=None)
#        data = net[net[0]!=net[1]]
#        data = data.reset_index()
        start = time.time()
        for i in range(len(data)):
            count += 1
            print("\r{} done".format(count), end="")
            #Alternative hypo
            x = data[0][i]
            y = data[1][i]
            if x==y:
                ofile.write("{}\t{}\t{}\n".format(x, y, data[2][i]))
                continue 
#            start = time.time()
            value = mast_dict[(x,y)]
        
            currentDf = pd.DataFrame.from_records([(x,y,cond,folder, value[cond][folder]['count'],
                                                value[cond][folder]['total'],value[cond][folder]['prob'],
                                                value[cond][folder]['pre']) for cond in value for folder in value[cond]])
            currentDf['key'] = 0
            matched = 0
            mergeDf = None
            for cond in conditions:
                if mergeDf is None:
                    mergeDf = currentDf[currentDf[2]==cond]
                else:
                    mergeDf = pd.merge(mergeDf,currentDf[currentDf[2]==cond], on='key')
            if  mergeDf is not None and len(mergeDf)>0:
                mergeDf['diff_prob'] = abs(mergeDf['6_x']-mergeDf['6_y'])
                mergeDf['threshold'] = mergeDf[['6_x','6_y']].max(axis=1)*0.8
                count1 = len(mergeDf[mergeDf['threshold']-mergeDf['diff_prob']>=0])
                if count1/7>=0.5: 
                    ofile.write("{}\t{}\t{}\n".format(x, y, data[2][i]))
#                    big_list.append(mergeDf.values.tolist())
#                if len(big_list) > 1000:
#                    print("Writing to file")
#                    write_to_file(big_list, big_file)                    
#wr.writerows(big_list)
                    big_list = []
     #               big_df = pd.concat([big_df, mergeDf])
        end = time.time()
        print("Total", end-start)
    
#    write_to_file(big_list, big_file)
    big_file.close()
#    big_df.to_csv(big_file, header=False)
    print("\nTrimmed {} edges".format(numTrimmed))

