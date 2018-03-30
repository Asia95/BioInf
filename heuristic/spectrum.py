import math
import random
import copy


temp = 30
temp2 = 32
graph = []
n = 0
sequence = ""


class Oligonucleotide:
    def __init__(self, i, series, size, min, max, first):
        self.id = i
        self.series = series
        self.size = size
        self.min = min
        self.max = max
        self.times_used = 0
        self.first = first
        self.komplementarny = False

    def __str__(self):
        return str(self.id) + ' : ' + self.series + ', ' + str(self.size) + ', ' + str(self.min) + ', ' + str(self.max) + ', ' + str(self.first)


class World:
    def __init__(self, oli):
        self.start = oli
        self.edges = []


class Edge:
    def __init__(self, end=None, pheromone=None):
        # self.start = None if start is None else start
        self.end = None if end is None else end
        self.move = []
        self.rest = []
        self.available = False
        self.pheromone = 0.1 if pheromone is None else pheromone

    def __str__(self):
        return str(self.rest) + ', ' + str(self.move) + ', ' + str(self.available)


def count_temp(sign):
    return {
        'A': 2,
        'T': 2,
        'C': 4,
        'G': 4,
    }[sign]


def temp_series(series):
    temp = 0
    for i in range(len(series)):
        temp = temp + count_temp(series[i])
    return temp


def complementary_sign(sign):
    return {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }[sign]


def spectrum():
    global graph
    global temp
    global temp2
    global n
    num_error = 1
    repeated_in = 3
    plik = open('seq.txt')
    try:
        f = plik.readlines()
        first_seq = f[0][:-1]
        #first_seq = first_seq[:-1]
        global sequence
        sequence = f[1][:-1]
        n = len(sequence)
    finally:
        plik.close()

    print(sequence)
    comp_sequence = ""
    for i in reversed(range(n)):
        comp_sequence = comp_sequence + complementary_sign(sequence[i])
    print(comp_sequence)

    elem_spectrum = []
    w = 0
    i = 0
    while i <= n:

        elem = sequence[w:w+i-w]
        if temp_series(elem) == temp or temp_series(elem) == temp2:
            elem_spectrum.append(elem)
            elem = comp_sequence[(len(comp_sequence)-i):(len(comp_sequence)-i)+(len(comp_sequence)-w-(len(comp_sequence)-i))]
            elem_spectrum.append(elem)
        i = i + 1
        if temp_series(elem) >= temp2:
            w = w + 1
            i = i - 1

    for p in elem_spectrum: print(p)

    repeated_times = 0
    tmp_id = 0
    for i in range(len(elem_spectrum)):
        repeated = False
        for j in range(len(graph)):
            if elem_spectrum[i] == graph[j].start.series:
                graph[j].start.min += 1
                repeated = True
                repeated_times += 1
        if repeated == False:
            if len(graph) == 0:
                tmp_id = 0
            else:
                tmp_id = graph[j].start.id + 1
            if elem_spectrum[i] == first_seq.strip():
                f = 1
            else:
                f = 0
            w = World(Oligonucleotide(tmp_id, elem_spectrum[i], len(elem_spectrum[i]), 1, 0, f))
            graph.append(w)
    hej = True
    #if (repeated_times/len(elem_spectrum)*100) >= repeated_in-1.5 and (repeated_times/len(elem_spectrum)*100) <= repeated_in+1.5:
    if hej == True:
        num_error2 = math.floor(len(elem_spectrum)*num_error/100)
        global missing
        missing = num_error2
        for i in range(num_error2):
            ok = False
            while ok == False:
                choosen = 0
                choosen = random.randint(0,len(graph))
                if graph[choosen].start.min != 0 and choosen != 0:
                    graph[choosen].start.min -= 1
                    ok = True

        i = 0
        while i < len(graph):

            graph[i].start.id = i

            if graph[i].start.min == 1:
                graph[i].start.min = 1
                graph[i].start.max = 1
            elif graph[i].start.min == 2 or graph[i].start.min == 3:
                graph[i].start.min = 2
                graph[i].start.max = 3
            elif graph[i].start.min == 4 or graph[i].start.min == 5:
                graph[i].start.min = 4
                graph[i].start.max = 5
            elif graph[i].start.min == 0:
                del graph[i]
                i = i - 1
            else:
                graph[i].start.min = -1
                graph[i].start.max = -1

            i = i + 1

        return graph
    else:
        return False

    for p in graph: print(p.start)


def return_len_seq():
    return n


def return_seq():
    return sequence

def matrix():
    for o in range(len(graph)):
        for g in range(len(graph)):
            count = 0
            vert = Edge()
            if len(graph[o].start.series) <= len(graph[g].start.series):
                vert.end = graph[g].start
                while count <= len(graph[o].start.series):
                    tmp = graph[o].start.series[count:count + len(graph[o].start.series) - count]
                    if tmp == graph[g].start.series[0:len(tmp)]:
                        if o != g or count != 0:
                            vert.move.append(count)
                            vert.rest.append(0)
                    count += 1
            else:
                vert.end = graph[g].start
                while count <= len(graph[o].start.series):
                    if len(graph[o].start.series) - count >= len(graph[g].start.series):
                        tmp = graph[o].start.series[count:count + len(graph[g].start.series)]
                    else:
                        tmp = graph[o].start.series[count:count + len(graph[o].start.series) - count]
                    if tmp == graph[g].start.series[0:len(tmp)]:
                        vert.move.append(count)
                        vert.rest.append(len(graph[o].start.series) - count - len(tmp))
                    count += 1
            vert.available = True
            vert.pheromones = 0
            graph[o].edges.append(vert)

    m = len(graph) - 1
    while m >= 0:
        mm = len(graph[m].edges) - 1
        while mm >= 0:
            mmm = len(graph[m].edges[mm].move) -1
            while mmm >= 0:
                if len(graph[m].start.series) == graph[m].edges[mm].move[mmm]:
                    del graph[m].edges[mm].move[mmm]
                    del graph[m].edges[mm].rest[mmm]
                mmm -= 1
            mm -= 1
        m -= 1
    m = len(graph) - 1
    while m >= 0:
        mm = len(graph[m].edges) - 1
        while mm >= 0:
            if len(graph[m].edges[mm].move) == 0:
                del graph[m].edges[mm]
            mm -= 1
        m -= 1


def complementary():
    global graph
    global missing
    for o in range(len(graph)):
        if graph[o].start.komplementarny is False:
            for j in range(len(graph)):
                if graph[j].start.komplementarny is False and len(graph[o].start.series) == len(graph[j].start.series):
                    komp = True
                    for k in range(len(graph[o].start.series)):
                        if graph[o].start.series[k] == 'A':
                            if graph[j].start.series[len(graph[j].start.series) - 1 - k] != 'T':
                                komp = False
                        elif graph[o].start.series[k] == 'T':
                            if graph[j].start.series[len(graph[j].start.series) - 1 - k] != 'A':
                                komp = False
                        elif graph[o].start.series[k] == 'C':
                            if graph[j].start.series[len(graph[j].start.series) - 1 - k] != 'G':
                                komp = False
                        elif graph[o].start.series[k] == 'G':
                            if graph[j].start.series[len(graph[j].start.series) - 1 - k] != 'C':
                                komp = False
                    if komp is True:
                        graph[o].start.komplementarny = True
                        graph[j].start.komplementarny = True
                        if graph[o].start.min != graph[j].start.min:
                            if graph[o].start.min == -1 or graph[j].start.min == -1:
                                graph[o].start.min = -1
                                graph[o].start.max = -1
                                graph[j].start.min = -1
                                graph[j].start.max = -1
                            elif graph[o].start.min > graph[j].start.min:
                                graph[j].start.min = graph[o].start.min
                                graph[j].start.max = graph[o].start.max
                            else:
                                graph[o].start.min = graph[j].start.min
                                graph[o].start.max = graph[j].start.max
                        break
    count = 0
    for o in range(len(graph)):
        if graph[o].start.komplementarny is False:
            missing -= 1
            new_series = []
            for i in range(len(graph[o].start.series)):
                new_series.append(complementary_sign(graph[o].start.series[i]))
            oli = copy.copy(graph[o].start)
            new_series = new_series[::-1]
            oli.series = ''.join(new_series)
            oli.id = len(graph)
            count += 1
            oli.komplementarny = True
            graph[o].start.komplementarny = True
            graph.append(World(copy.copy(oli)))