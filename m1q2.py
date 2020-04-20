import networkx as nx
import EoN
import matplotlib.pyplot as plt
import numpy as np


N = 30
edges = 5
realisations = 10
tmax=20
tmin=0
G = nx.watts_strogatz_graph(N, edges, 0) #Watts-Strogatz Small World Graph

beta = [0.25, 0.50, 0.75]
gamma = 0.5


def properRound(num, dec=0):
    num = str(num)[:str(num).index('.')+dec+2]
    if num[-1]>='5':
      a = num[:-2-(not dec)]       # integer part
      b = int(num[-2-(not dec)])+1 # decimal part
      return float(a)+b**(-dec+1) if a and b == 10 else float(a+str(b))
    return float(num[:-1])

def sortArray(array):
    rightArray = array[len(array)//2: len(array)+1]
    for element in range(len(rightArray)):
        array.remove(rightArray[element])
        array.insert(element, rightArray[element])
    return array


def getAllData(G = G, N = N, beta = beta, gamma = gamma, realisations = realisations, tmin = tmin, tmax = tmax, extraEdge = False):
    if(extraEdge):
        G.add_edge(0, N//2)
    #nx.draw(G, pos=nx.circular_layout(G), node_color='r', edge_color='b')

    dAllData = {}
    for node in range(N):
        dBetas = {}
        for be in beta:
            dRealisations = {}
            for rs in range(realisations):
                t, S, I, R = EoN.fast_SIR(G, be, gamma, initial_infecteds = [node], tmin=tmin, tmax=tmax )
                dResults={}
                for index in range(len(t)):
                    time= (t[index])
                    if time in dResults:
                        dResults[time] = int(round(dResults.get(time) + I[index]) / 2)
                    else:
                        dResults[time] = int(I[index])
                dRealisations["R" + str(rs+1)] = dResults
            dBetas[str(be)] = dRealisations
        dAllData[str(node+1)] = dBetas
    return dAllData


def getAvgOfRealisations(data):
    for nodes in data:
        betas = data[nodes]
        for beta in betas:
            eachRealisations = {}
            realisations = betas[beta]
            for realisation in realisations:
                singleRealisation = realisations[realisation]
                for time in singleRealisation:
                    if time in eachRealisations:
                        eachRealisations[time] = (eachRealisations.get(time) + singleRealisation[time]) / 2
                    else:
                        eachRealisations[time] = singleRealisation[time]

            sortedEachRealisations = {}
            for key in sorted(eachRealisations):
                sortedEachRealisations[key] = eachRealisations.get(key)
            data[nodes][beta] = {sum(sortedEachRealisations.keys())/len(sortedEachRealisations): sum(sortedEachRealisations.values())/len(sortedEachRealisations)}
    return data


def getMaxOfRealisations(data):
    for nodes in data:
        betas = data[nodes]
        for beta in betas:
            eachRealisations = {}
            realisations = betas[beta]
            for realisation in realisations:
                singleRealisation = realisations[realisation]
                for time in singleRealisation:
                    if time in eachRealisations:
                        eachRealisations[time] = (eachRealisations.get(time) + singleRealisation[time]) / 2
                    else:
                        eachRealisations[time] = singleRealisation[time]

            sortedEachRealisations = {}
            for key in sorted(eachRealisations):
                sortedEachRealisations[key] = eachRealisations.get(key)
            maxKey = max(sortedEachRealisations, key=sortedEachRealisations.get)
            data[nodes][beta] = {maxKey: sortedEachRealisations.get(maxKey)}
    return data


def getAllBetas(data, beta):
    nodesWith25BetaTime = []
    nodesWith25BetaValues = []
    nodeName = []
    for nodes in data:
        nodeName.append(nodes)
        time = data[nodes][str(beta)]
        for values in time:
            nodesWith25BetaTime.append(values)
            nodesWith25BetaValues.append(time[values])

    return nodesWith25BetaTime, nodesWith25BetaValues, nodeName


def getBarPlot(data, dataExtraEdge, beta, sortAlgorithm):
    time, values, nodes = data
    timeExtraEdge, valuesExtraEdge, nodesExtraEdge = dataExtraEdge

    nodes = sortAlgorithm(nodes)
    values = sortAlgorithm(values)
    valuesExtraEdge = sortAlgorithm(valuesExtraEdge)

    x = np.arange(len(nodes))
    width = 0.35

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, values, width, label="Normal")
    rects2 = ax.bar(x + width/2, valuesExtraEdge, width, label="ExtraEdge")

    ax.set_ylabel("#InfectedPeople")
    ax.set_xlabel("Nodes")
    ax.set_title("Max number of Infected People per initial Node with beta = " + str(beta))
    ax.set_xticks(x)
    ax.set_xticklabels(nodes)
    ax.legend()

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')


    autolabel(rects1)
    autolabel(rects2)
    fig.tight_layout()
    plt.show()


def getScatterPlot(data, dataExtraEdge, beta):
    time, values, nodes = data
    timeExtraEdge, valuesExtraEdge, nodesExtraEdge = dataExtraEdge

    fig, ax = plt.subplots()
    ax.scatter(values, time, label="Normal")
    ax.scatter(valuesExtraEdge, timeExtraEdge, label="ExtraEdge")

    for i, nodeName in enumerate(nodes):
        ax.annotate(nodeName, (values[i], time[i]))

    for i, nodeName in enumerate(nodesExtraEdge):
        ax.annotate(nodeName, (valuesExtraEdge[i], timeExtraEdge[i]))

    plt.title("beta = " + str(beta))
    plt.ylabel("time")
    plt.xlabel("Infected People")
    plt.show()


def compareResultsInDifferentBetas(maxInfextedPeople25, maxInfextedPeople50, maxInfextedPeople75, sortAlgorithm, hasExtraEdge = False):
    message = " "
    if(hasExtraEdge):
        message = "with an extra edge"

    time25, values25, nodes25 = maxInfextedPeople25
    time50, values50, nodes50 = maxInfextedPeople50
    time75, values75, nodes75 = maxInfextedPeople75

    nodes25 = sortAlgorithm(nodes25)
    values25 = sortAlgorithm(values25)
    values50 = sortAlgorithm(values50)
    values75 = sortAlgorithm(values75)

    x = np.arange(len(nodes25))
    width = 0.30

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, values25, width, label="Beta = 25 " + message)
    rects2 = ax.bar(x + width/2, values50, width, label="Beta = 50" + message)
    rects3 = ax.bar(x + 3*(width/2), values75, width, label="Beta = 75" + message)

    ax.set_ylabel("#InfectedPeople")
    ax.set_xlabel("Nodes")
    ax.set_title("Max number of Infected People per initial Node with different beta")
    ax.set_xticks(x)
    ax.set_xticklabels(nodes25)
    ax.legend()

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    fig.tight_layout()
    plt.show()


#compareResultsInDifferentBetas(getAllBetas(getMaxOfRealisations(getAllData()),  0.25), getAllBetas(getMaxOfRealisations(getAllData()),  0.5),
#    getAllBetas(getMaxOfRealisations(getAllData()),  0.75), sortArray)

#compareResultsInDifferentBetas(getAllBetas(getMaxOfRealisations(getAllData(extraEdge=True)),  0.25), getAllBetas(getMaxOfRealisations(getAllData(extraEdge=True)),  0.5),
#    getAllBetas(getMaxOfRealisations(getAllData(extraEdge=True)),  0.75), sortArray, hasExtraEdge=True)

getScatterPlot(getAllBetas(getMaxOfRealisations(getAllData()),  0.75), getAllBetas(getMaxOfRealisations(getAllData(extraEdge=True)),  0.75), 0.75 )
#getBarPlot(getAllBetas(getMaxOfRealisations(getAllData()),  0.25), getAllBetas(getMaxOfRealisations(getAllData(extraEdge=True)),  0.25), 0.25, sortArray )
