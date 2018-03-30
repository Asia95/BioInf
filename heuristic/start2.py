import math
import random
import time
import copy
from spectrum import *


class Pheromon:
    def __init__(self):
        self.start = None
        self.end = None
        self.phero = 0


class Seq:
    def __init__(self):
        self.chain = ""
        self.correct = False


def get_complementary(spectrum, id):
    tmp = ""
    for o in range(len(spectrum[id].start.series)):
        tmp += (complementary_sign(spectrum[id].start.series[o]))

    tmp = tmp[::-1]
    for o in range(len(spectrum)):
        if tmp == spectrum[o].start.series:
            return o
    return -1


def set_available(id1, id2, graph):
    if graph[id1].start.max != -1:
        if graph[id1].start.times_used >= graph[id1].start.max:
            for s in range(len(graph)):
                for ss in range(len(graph[s].edges)):
                    if graph[s].edges[ss].end.id == id1:
                        graph[s].edges[ss].available = False
                    if graph[s].edges[ss].end.id == id2:
                        graph[s].edges[ss].available = False
        else:
            for s in range(len(graph)):
                for ss in range(len(graph[s].edges)):
                    if graph[s].edges[ss].end.id == id1:
                        graph[s].edges[ss].available = True
                    if graph[s].edges[ss].end.id == id2:
                        graph[s].edges[ss].available = True
    return graph


def add_olig(i, count, spectrum, sequence2):
    tmp = spectrum[i].start.series[0:count]
    for j in range(len(tmp)):
        sequence2.append(tmp[j])
    return sequence2


def add_series_end(i, j, spectrum, series):
    tmp = spectrum[j].start.series
    for t in range(len(tmp)):
        series.insert(i, tmp[t])
        i += 1
    return series


def quality(seq):
    global sequence
    good = 0
    for i in range(len(seq)):
        if sequence[i] == seq[i]:
            good += 1
    return good / len(seq) * 100


count2 = 0


def check_if_ok(spectrum, seq):
    occur = []
    found = False
    global count2

    for o in range(len(spectrum)):
        if spectrum[o].start.times_used < spectrum[o].start.min: # spectrum[o].start.times_used > spectrum[o].start.max:
            return False

    new_seq = ''.join(seq)

    if len(seq) > n:
        new_seq_correct = False
    else:
        new_seq_correct = True
    """
    for s in range(len(graph)):
        occur.clear()
        id1 = new_seq.find(graph[s].start.series)
        if id1 != -1:
            occur.append(id1)
        while id1 >= 0:
            id1 = new_seq.find(graph[s].start.series, id1 + 1)
            if id1 != -1:
                occur.append(id1)

        if graph[s].start.min == 0:
            if len(occur) > 0 and len(occur) < 1:
                new_seq_correct = False
        elif graph[s].start.min == 1:
            if len(occur) > 1 and len(occur) < 3:
                new_seq_correct = False
        elif graph[s].start.min == 2:
            if len(occur) > 2 and len(occur) < 5:
                new_seq_correct = False
        elif graph[s].start.min == 4:
            if len(occur) > 4:
                new_seq_correct = False
        elif graph[s].start.min == -1:
            if len(occur) >= 5:
                new_seq_correct = False
    """

    if new_seq_correct == True:
        count2 += 1
        print(str(count2))
        print("wygenerowana sekwencja: " + new_seq)
        print("jakość dopasowania: " + '%.2f' % quality(new_seq) + " %")
        print()
        found = True
    if found == True:
        return True
    else:
        return False


def ant(graph, first):
    firstCount = True
    ants = 20
    iterac = 60
    par = 20
    divide = 90
    pheromoneFactor = 30
    minpheromone = 10
    choosenRoads = []
    pheroToAdd = []
    okSeq = []
    sequence2 = []
    print(str(n))
    for i in range(iterac):
        for j in range(ants):
            # use first
            actual = first
            graph[actual].start.times_used += 1
            idx = get_complementary(graph, actual)
            graph[idx].start.times_used += 1
            graph = set_available(idx, actual, graph)
            # add those who are in first
            for k in range(len(graph[actual].edges)):
                for m in range(len(graph[actual].edges[k].rest)):
                    if graph[actual].edges[k].rest[m] > 0:
                        graph[graph[actual].edges[k].end.id].start.times_used += 1
                        idx = get_complementary(graph, k)
                        graph[idx].start.times_used += 1
                        graph = set_available(idx, k, graph)
            minOli = False
            roads = True
            tooLong = False
            while minOli == False and roads == True and tooLong == False:
                if firstCount:
                    # choose road

                    for k in range(len(graph[actual].edges)):
                        # print(str(k) + " " + str(len(graph[actual].edges)))
                        if graph[graph[actual].edges[k].end.id].start.times_used < graph[graph[actual].edges[k].end.id].start.max:
                            if graph[actual].edges[k].rest[0] == 0:
                                choosenRoads.append(graph[actual].edges[k].end.id)
                                if graph[graph[actual].edges[k].end.id].start.times_used < graph[graph[actual].edges[k].end.id].start.min:
                                    choosenRoads.append(graph[actual].edges[k].end.id)
                                    choosenRoads.append(graph[actual].edges[k].end.id)
                else :
                    # choose road with pheromones
                    for k in range(len(graph[actual].edges)):
                        if graph[graph[actual].edges[k].end.id].start.times_used < graph[graph[actual].edges[k].end.id].start.max and graph[actual].edges[k].rest[0] == 0:
                            choosenRoads.append(graph[actual].edges[k].end.id)
                            for m in range(graph[actual].edges[k].pheromones * pheromoneFactor):
                                choosenRoads.append(graph[actual].edges[k].end.id)
                            if graph[graph[actual].edges[k].end.id].start.times_used < graph[graph[actual].edges[k].end.id].start.min:
                                for m in range(minpheromone):
                                    choosenRoads.append(graph[actual].edges[k].end.id)
                if len(choosenRoads) != 0:
                    # random road
                    ran = random.randint(0,len(choosenRoads))
                    # print("ran : " + str(ran) + " , " + str(len(choosenRoads)))
                    if ran == len(choosenRoads) and ran != 0:
                        # print("ran : " + str(ran) + " , " + str(len(choosenRoads)))
                        ran -= 1
                    next = choosenRoads[ran]
                    found = False
                    find_next = 0
                    for hej in range(len(graph[actual].edges)):
                        if next == graph[actual].edges[hej].end.id and found == False:
                            found = True
                            find_next = hej
                            break
                            #print(str(hej))
                    """
                    for index, item in enumerate(graph[actual].edges):
                        if next == item.end.id:
                            next2 = index
                            #break
                    """
                    graph[next].start.times_used += 1
                    idx = get_complementary(graph, next)
                    graph[idx].start.times_used += 1
                    graph = set_available(idx, next, graph)
                    # add those who are in next
                    for k in range(len(graph[next].edges)):
                        for m in range(len(graph[next].edges[k].rest)):
                            if graph[next].edges[k].rest[m] > 0:
                                graph[graph[next].edges[k].end.id].start.times_used += 1
                                idx = get_complementary(graph, graph[next].edges[k].end.id)
                                graph[idx].start.times_used += 1
                                graph = set_available(idx, graph[next].edges[k].end.id, graph)
                    # add to seq

                    sequence2 = add_olig(actual, graph[actual].edges[find_next].move[0], graph, sequence2)
                    if len(sequence2) > n:
                        tooLong = True
                    else:
                        # remeber and leave pheromone
                        pheromone = Pheromon()
                        pheromone.start = actual
                        pheromone.end = next
                        pheromone.phero = math.floor(divide / (graph[actual].edges[find_next].move[0] + 1))
                        pheroToAdd.append(pheromone)
                        # new vertex
                        actual = next
                        choosenRoads = []
                        # chceck if minOli is used
                        for k in range(len(graph)):
                            if graph[k].start.times_used < graph[k].start.min:
                                minOli = False
                            else:
                                minOli = True
                else:
                    roads = False
            # add end
            sequence2 = add_series_end(len(sequence2), actual, graph, sequence2)
            # check if sequence is ok
            if len(sequence2) <= n:
                if check_if_ok(graph, sequence2) == True:
                    okSeq.append(sequence2)
                #else:
                    #print("wrong: " + ''.join(sequence2))
            #else:
                #print("wrong2: " + ''.join(sequence2))
            sequence2 = []
            for k in range(len(graph)):
                graph[k].start.times_used = 0
        # add pheromone
        for j in range(len(pheroToAdd)):
            for jj in range(len(graph[pheroToAdd[j].start].edges)):
                if graph[pheroToAdd[j].start].edges[jj].end.id == pheroToAdd[j].end:
                    graph[pheroToAdd[j].start].edges[jj].pheromones = pheroToAdd[j].phero
                    break
        for j in range(len(graph)):
            for k in range(len(graph[j].edges)):
                if graph[j].edges[k].pheromones > 0:
                    if graph[j].edges[k].pheromones < par:
                        graph[j].edges[k].pheromones = 0
                    else:
                        graph[j].edges[k].pheromones = graph[j].edges[k].pheromones - par
        pheroToAdd = []
        firstCount = False
    # for i in range(len(okSeq)):
        # print("Sequence: ")
        # print(okSeq[i])


if __name__ == '__main__':

    graph = spectrum()
    n = return_len_seq()
    sequence = return_seq()
    complementary()
    matrix()

    start = time.clock()

    ant(graph, 0)

    print("time: " + str(time.clock() - start))
    print("konice")