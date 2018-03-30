import math
import random

temp = 30
temp2 = 32
n = 0
sequence = ""
graph = []


class Oligonucleotide:
    def __init__(self, i, series, size, min, max, first):
        self.id = i
        self.series = series
        self.size = size
        self.min = min
        self.max = max
        self.max2 = 0
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
    }.get(sign, 'X')


def spectrum():
    global graph
    global temp
    global temp2
    global n
    num_error = 1
    # num_repeated = 1
    plik = open('seq.txt')
    try:
        f = plik.readlines()
        first_seq = f[0]
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
    # print(len(comp_sequence))

    elem_spectrum = []
    w = 0
    i = 0
    while i <= n:

        elem = sequence[w:w+i-w]
        # print("elem: " + elem)
        # print("temp: " + str(temp_series(elem)))
        if temp_series(elem) == temp or temp_series(elem) == temp2:
            elem_spectrum.append(elem)
            elem = comp_sequence[(len(comp_sequence)-i):(len(comp_sequence)-i)+(len(comp_sequence)-w-(len(comp_sequence)-i))]
            # print("elem2: " + elem)
            elem_spectrum.append(elem)
        i = i + 1
        if temp_series(elem) >= temp2:
            w = w + 1
            i = i - 1

    for p in elem_spectrum: print(p)

    # print("elem_spectrum : " + str(len(elem_spectrum)))
    # print("spectrum2 : " + str(len(spectrum2)))
    repeated_times = 0
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
                # tmp_id = len(spectrum2) + 1
                tmp_id = graph[j].start.id + 1
            if elem_spectrum[i] == first_seq.strip():
                f = 1
            else:
                f = 0
            w = World(Oligonucleotide(tmp_id, elem_spectrum[i], len(elem_spectrum[i]), 1, 0, f))
            graph.append(w)

    # print(str(repeated_times))

    if (repeated_times/len(elem_spectrum)*100) >= -0.5 and (repeated_times/len(elem_spectrum)*100) <= 2.5:

        # print("ok")
        # ok = False
        # choosen = 0
        num_error2 = math.floor(len(elem_spectrum)*num_error/100)
        global missing
        missing = num_error2
        for i in range(num_error2):
            ok = False
            while ok == False:
                choosen = random.radiant(0,len(graph))
                # print(str(choosen))
                if graph[choosen].start.min != 0:
                    # print(str(choosen))
                    graph[choosen].start.min -= 1
                    ok = True

        i = 0
        # print(first_seq)
        while i < len(graph):

            """if spectrum2[i].min == 1:
                spectrum2[i].min = 1
                spectrum2[i].max = 1
                spectrum2[i].max2 = 3
            elif spectrum2[i].min == 2 or spectrum2[i].min == 3:
                spectrum2[i].min = 2
                spectrum2[i].max = 3
                spectrum2[i].max2 = 5
            elif spectrum2[i].min == 4 or spectrum2[i].min == 5:
                spectrum2[i].min = 4
                spectrum2[i].max = 5
                spectrum2[i].max2 = -1
            elif spectrum2[i].min == 0:
                del spectrum2[i]
                i = i - 1
            else:
                spectrum2[i].min = -1
                spectrum2[i].max = -1
                spectrum2[i].max2 = -1"""

            if graph[i].start.min == 1 or graph[i].start.min == 2:
                graph[i].start.min = 1
                graph[i].start.max = 2
                graph[i].start.max2 = 4
            elif graph[i].start.min == 3 or graph[i].start.min == 4:
                graph[i].start.min = 3
                graph[i].start.max = 4
                graph[i].start.max2 = -1
            elif graph[i].start.min == 0:
                del graph[i]
                i = i - 1
            else:
                graph[i].start.min = -1
                graph[i].start.max = -1
                graph[i].start.max2 = -1

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
    # print(str(len(spectrum)))
    for o in range(len(graph)):
        for g in range(len(graph)):
            count = 0
            # ok = False
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
            graph[o].edges.append(vert)
            # vert.move[:] = []
            # vert.rest = []

    ile = 0
    for o in range(len(graph)):
        if temp_series(graph[o].start.series) == temp:
            ile = ile + graph[o].start.min
    # print(str(ile))
    # print("matrix move: " + str((spec_matrix[0][1].move)))
    m = len(graph) - 1
    while m >= 0:
        mm = len(graph[m].edges) - 1
        while mm >= 0:
            mmm = len(graph[m].edges[mm].move) - 1
            while mmm >= 0:
                if (ile/2) - 2 + graph[m].edges[mm].move[mmm] + len(graph[m].start.series) > n and len(graph[m].edges[mm].move) > 0:
                    del graph[m].edges[mm].move[mmm:len(graph[m].edges[mm].move)]
                    del graph[m].edges[mm].rest[mmm:len(graph[m].edges[mm].rest)]
                    break
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
                                graph[o].start.max2 = -1
                                graph[j].start.min = -1
                                graph[j].start.max = -1
                                graph[j].start.max2 = -1
                            elif graph[o].start.min > graph[j].start.min:
                                graph[j].start.min = graph[o].start.min
                                graph[j].start.max = graph[o].start.max
                                graph[j].start.max2 = graph[o].start.max2
                            else:
                                graph[o].start.min = graph[j].start.min
                                graph[o].start.max = graph[j].start.max
                                graph[o].start.max2 = graph[j].start.max2
                        break
    count = 0
    for o in range(len(graph)):
        if graph[o].start.komplementarny is False:
            missing -= 1
            new_series = ""
            for i in range(len(graph)):
                new_series.append(complementary_sign(graph[o].start.series[i]))
            oli = graph[o].start
            new_series = new_series[::-1]
            oli.series = new_series
            oli.id = len(graph) + count
            count += 1
            oli.komplementarny = True
            graph[o].start.komplementarny = True
            graph.append(World(list(oli)))