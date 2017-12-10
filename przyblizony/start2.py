import math
import random
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


def add_olig(i, count, spectrum, sequence):
    tmp = spectrum[i].start.series[0:count]
    for j in range(len(tmp)):
        sequence.append(tmp[j])
    return sequence


def add_series_end(i, j, spectrum, series):
    tmp = spectrum[j].start.series
    for t in range(len(tmp)):
        series.insert(i, tmp[t])
        i += 1
    return series


def quality(seq):
    global sequence
    good = 0
    for i in range(n):
        if sequence[i] == seq[i]:
            good += 1
    return good / n * 100


count2 = 0


def check_if_ok(spectrum, seq):
    occur = []
    found = False
    global count2

    for o in range(len(spectrum)):
        if spectrum[o].start.times_used < spectrum[o].start.min and spectrum[o].start.times_used > spectrum[o].start.max:
            return False

    new_seq = ''.join(seq)

    if len(new_seq) != n:
        new_seq_correct = False
    else:
        new_seq_correct = True

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
            if len(occur) > 0 and len(occur) < 2:
                new_seq_correct = False
        elif graph[s].start.min == 1:
            if len(occur) > 1 and len(occur) < 4:
                new_seq_correct = False
        elif graph[s].start.min == 3:
            if len(occur) > 3:
                new_seq_correct = False
        elif graph[s].start.min == -1:
            if len(occur) >= 4:
                new_seq_correct = False

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
    ants = 50
    iterac = 30
    par = 20
    divide = 90
    pheromoneFactor = 20
    minpheromone = 10
    choosenRoads = []
    pheroToAdd = []
    okSeq = []
    sequence = []

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
                    for k in range(len(graph[actual].edges)):
                        if next == graph[actual].edges[k].end.id:
                            next2 = k
                            break
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
                    sequence = add_olig(actual, graph[actual].edges[next2].move[0], graph, sequence)
                    if len(sequence) > n:
                        tooLong = True
                    else:
                        # remeber and leave pheromone
                        pheromone = Pheromon()
                        pheromone.start = actual
                        pheromone.end = next
                        pheromone.phero = math.floor(divide / (graph[actual].edges[next2].move[0] + 1))
                        pheroToAdd.append(pheromone)
                        # new vertex
                        actual = next
                        choosenRoads = []
                        # chceck if minOli is used
                        for k in range(len(graph)):
                            if graph[k].start.times_used < graph[k].start.min:
                                minOli = False
                                break
                            minOli = True
                else:
                    roads = False
            # add end
            sequence = add_series_end(len(sequence), actual, graph, sequence)
            # check if sequence is ok
            if len(sequence) <= n:
                if check_if_ok(graph, sequence) == True:
                    okSeq.append(sequence)
                else:
                    print("erong: " + sequence)
            sequence = []
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
    for i in range(len(okSeq)):
        print("Sequence: ")
        print(okSeq[i])


if __name__ == '__main__':

    spectrum()
    n = return_len_seq()
    sequence = return_seq()
    complementary()
    matrix()
    ant(graph, 0)
    print("konice")