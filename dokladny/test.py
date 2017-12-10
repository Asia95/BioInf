import string
import random
import math

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

    def get_series(self):
        return self.series

    def get_id(self):
        return self.id

    def get_min(self):
        return self.min

    def __str__(self):
        return str(self.id) + ' : ' + self.series + ', ' + str(self.size) + ', ' + str(self.min) + ', ' + str(self.max) + ', ' + str(self.first)
"""
class Vertex:
    def __init__(self, node):
        self.id = node
        self.neighbors = {}

    def __str__(self):
        return str(self.id.get_series()) + ' : ' + str([x.id.get_series()  + ' : ' + str(self.get_weight(x)) for x in self.neighbors])

    def __iter__(self):
        return iter(self.neighbors.values())

    def add_neighbor(self, neighbor, weight):
        self.neighbors[neighbor] = weight

    def get_neighbors(self):
        return self.neighbors.keys()

    def get_id(self):
        return self.id

    def get_weight(self, neighbor):
        return self.neighbors[neighbor]

    #czy końcówka taka sama jak początek kolejnego
    def endBeginningMatch(self, right):
        t = self.get_id().get_series()
        r = right.get_id().get_series()
        print(t)
        for i in range(len(t)):
            if t[i:] == r[:len(t)-i]:
                return t[i:]
        return ''

    def find_weight(self, next):
        tmp = 0
        if len(self.endBeginningMatch(next)) > 0:
            tmp = len(self.endBeginningMatch(next))
        else:
            if self.id.get_series() in next:
                tmp = len(self.id.get_series())
        return tmp

    def findNeighbourWithLargestWeight(self):
        m = 0
        for n in self.neighbors:
            if self.get_weight(n) > m:
                m = self.get_weight(n)
                v = n
        return v

class Graph:
    def __init__(self):
        self.vert_dict = {}
        self.num_vertices = 0

    def __iter__(self):
        return iter(self.vert_dict.values())

    def add_vertex(self, node):
        self.num_vertices = self.num_vertices + 1
        new_vertex = Vertex(node)
        self.vert_dict[node] = new_vertex
        return new_vertex

    def get_vertex(self, n):
        if n in self.vert_dict:
            return self.vert_dict[n]
        else:
            return None

    def add_edge(self, frm, to, cost):
        if frm not in self.vert_dict:
            self.add_vertex(frm)
        if to not in self.vert_dict:
            self.add_vertex(to)

        self.vert_dict[frm].add_neighbor(self.vert_dict[to], cost)
        #self.vert_dict[to].add_neighbor(self.vert_dict[frm], cost)

    def get_vertices(self):
        return self.vert_dict.keys()

def printInMatrix(g):
    print ('        ',end=" ")
    for j in sorted(g, key=id):
        print ("% 3s" % j.id,end=" ")
    print ()

    for i in sorted(g, key=id):
        print ("% 3s" % i.id,end=" ")
        for j in sorted(g, key=id):
            if i.get_id == j.get_id:
                s = '  -'
            else:
                s = "% 3s" % i.get_weight(j)
            print (s,end="      ")
        print ()
"""
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
"""
def first_spectrum(spectrum2):
    plik = open('plik.txt')
    try:
        f = plik.readlines()
        i = 0
        for line in f:
            word = line.split(' ')
            o = Oligonucleotide(i, word[0], len(word[0]), word[1], word[2], word[3])
            i = i + 1
            spectrum2.append(o)
            print(o)
    finally:
        plik.close()
"""
def spectrum(spectrum2):
    tmp = 30
    tmp2 = 32
    num_error = 1
    num_repeated = 1
    plik = open('seq.txt')
    try:
        f = plik.readlines()
        first_seq = f[0];
        sequence = f[1]
        n = len(sequence)-1
        #print(n)
    finally:
        plik.close()

    print(sequence)
    comp_sequence = ""
    for i in reversed(range(n)):
        comp_sequence = comp_sequence + complementary_sign(sequence[i])
    print(comp_sequence)
    #print(len(comp_sequence))

    elem_spectrum = []
    w = 0
    i = 0
    while i <= n:

        elem = sequence[w:w+i-w]
        #print("elem: " + elem)
        #print("temp: " + str(temp_series(elem)))
        if temp_series(elem) == tmp or temp_series(elem) == tmp2:
            elem_spectrum.append(elem)
            elem = comp_sequence[(len(comp_sequence)-i):(len(comp_sequence)-i)+(len(comp_sequence)-w-(len(comp_sequence)-i))]
            #print("elem2: " + elem)
            elem_spectrum.append(elem)
        i = i + 1
        if temp_series(elem) >= tmp2:
            w = w + 1
            i = i - 1

    for p in elem_spectrum: print(p)

    #print("elem_spectrum : " + str(len(elem_spectrum)))
    #print("spectrum2 : " + str(len(spectrum2)))
    repeated_times = 0
    for i in range(len(elem_spectrum)):
        repeated = False
        for j in range(len(spectrum2)):
            if elem_spectrum[i] == spectrum2[j].get_series():
                spectrum2[j].min += 1
                repeated = true
                repeated_times += 1
        if repeated == False:
            if len(spectrum2) == 0:
                tmp_id = 0
            else:
                tmp_id = spectrum2[j].get_id() + 1
            if elem_spectrum[i] == first_seq.strip():
                f = 1
            else:
                f = 0
            o = Oligonucleotide(tmp_id, elem_spectrum[i], len(elem_spectrum[i]), 1, 0, f)
            spectrum2.append(o)

    #print(str(repeated_times))

    if (repeated_times/len(elem_spectrum)*100) >= -0.5 and (repeated_times/len(elem_spectrum)*100) <= 2.5:

        #print("ok")
        ok = False
        choosen = 0
        num_error2 = math.floor(len(elem_spectrum)*num_error/100)
        for i in range(num_error2):
            ok = False
            while ok == False:
                choosen = random.radiant(0,len(spectrum2))
                #print(str(choosen))
                if spectrum2[choosen].min != 0:
                    #print(str(choosen))
                    spectrum2[choosen].min -= 1
                    ok = true

        i = 0
        #print(first_seq)
        while i < len(spectrum2):

            if spectrum2[i].min == 1:
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
                spectrum2[i].max2 = -1

            i = i + 1

    for p in spectrum2: print(p)

if __name__ == '__main__':

    spectrum2 = []

    #first_spectrum(spectrum2)
    spectrum(spectrum2)
    """
    g = Graph()

    for i in range(len(spectrum2)):
        g.add_vertex(spectrum2[i])

    for v in g:
        for w in g:
            if v.get_id() is not w.get_id():
                #g.add_edge(v.get_id(), w.get_id(), len(v.endBeginningMatch(w.get_id())))
                g.add_edge(v.get_id(), w.get_id(), v.find_weight(w.get_id()))

    for v in g:
        print (g.vert_dict[v.get_id()])

    print()
    printInMatrix(g)
"""
