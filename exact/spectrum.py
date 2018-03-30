import math
import random
import copy

temp = 30
temp2 = 32
n = 0
sequence = ""
spectrum2 = []


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


class Vertex:
    def __init__(self):
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
    global spectrum2
    tmp = 30
    tmp2 = 32
    num_error = 1
    repeated_in = 1
    # num_repeated = 1
    plik = open('seq.txt')
    try:
        f = plik.readlines()
        first_seq = f[0]
        global sequence
        sequence = f[1][:-1]
        global n
        n = len(sequence)
        # print(n)
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
        if temp_series(elem) == tmp or temp_series(elem) == tmp2:
            elem_spectrum.append(elem)
            elem = comp_sequence[(len(comp_sequence)-i):(len(comp_sequence)-i)+(len(comp_sequence)-w-(len(comp_sequence)-i))]
            # print("elem2: " + elem)
            elem_spectrum.append(elem)
        i = i + 1
        if temp_series(elem) >= tmp2:
            w = w + 1
            i = i - 1

    for p in elem_spectrum: print(p)

    # print("elem_spectrum : " + str(len(elem_spectrum)))
    # print("spectrum2 : " + str(len(spectrum2)))
    repeated_times = 0
    for i in range(len(elem_spectrum)):
        repeated = False
        for j in range(len(spectrum2)):
            if elem_spectrum[i] == spectrum2[j].series:
                spectrum2[j].min += 1
                repeated = True
                repeated_times += 1
        if repeated == False:
            if len(spectrum2) == 0:
                tmp_id = 0
            else:
                # tmp_id = len(spectrum2) + 1
                tmp_id = spectrum2[j].id + 1
            if elem_spectrum[i] == first_seq.strip():
                f = 1
            else:
                f = 0
            o = Oligonucleotide(tmp_id, elem_spectrum[i], len(elem_spectrum[i]), 1, 0, f)
            spectrum2.append(o)

    print(str(repeated_times))
    print(str(repeated_times/len(elem_spectrum)*100))

    if (repeated_times/len(elem_spectrum)*100) >= repeated_in - 1.5 and (repeated_times/len(elem_spectrum)*100) <= repeated_in + 1.5:

        # print("ok")
        # ok = False
        # choosen = 0
        num_error2 = math.floor(len(elem_spectrum)*num_error/100)
        global missing
        missing = num_error2
        for i in range(num_error2):
            ok = False
            while ok == False:
                choosen = 0
                choosen = random.randint(0,len(spectrum2))
                # print(str(choosen))
                if spectrum2[choosen].min != 0:
                    # print(str(choosen))
                    spectrum2[choosen].min -= 1
                    ok = True

        i = 0
        # print(first_seq)
        while i < len(spectrum2):

            spectrum2[i].size = len(spectrum2[i].series)
            spectrum2[i].times_used = 0
            spectrum2[i].komplementarny = False
            spectrum2[i].id = i

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

        return spectrum2
    else:
        return False


def return_len_seq():
    return n


def return_seq():
    return sequence


def matrix(spectrum):
    spec_matrix = []
    # print(str(len(spectrum)))
    for o in range(len(spectrum)):
        col = []
        for g in range(len(spectrum)):
            count = 0
            # ok = False
            vert = Vertex()
            if len(spectrum[o].series) <= len(spectrum[g].series):
                while count <= len(spectrum[o].series):
                    tmp = spectrum[o].series[count:count + len(spectrum[o].series) - count]
                    if tmp == spectrum[g].series[0:len(tmp)]:
                        if o != g or count != 0:
                            vert.move.append(count)
                            vert.rest.append(0)
                    count += 1
            else:
                while count <= len(spectrum[o].series):
                    if len(spectrum[o].series) - count >= len(spectrum[g].series):
                        tmp = spectrum[o].series[count:count + len(spectrum[g].series)]
                    else:
                        tmp = spectrum[o].series[count:count + len(spectrum[o].series) - count]
                    if tmp == spectrum[g].series[0:len(tmp)]:
                        vert.move.append(count)
                        vert.rest.append(len(spectrum[o].series) - count - len(tmp))
                    count += 1
            vert.available = True
            col.append(vert)
            # vert.move[:] = []
            # vert.rest = []
        # print(str(col[0][0].move[0]))
        spec_matrix.append(list(col))
        # col.clear()

    ile = 0
    # ile2 = 0
    # n = 20
    for o in range(len(spectrum)):
        if temp_series(spectrum[o].series) == temp:
            ile = ile + spectrum[o].min
    # print(str(ile))
    # print("matrix move: " + str((spec_matrix[0][1].move)))

    for m in range(len(spec_matrix)):
        for mm in range(len(spec_matrix[m])):
            # print("matrix move: " + str(len(spec_matrix[m][mm].move)))
            for k in range(len(spec_matrix[m][mm].move)):
                # print("m: " + str(m) + " mm: " + str(mm) + " k: " + str(k))
                # print("matrix move: " + str(len(spec_matrix[m][mm].move)))
                # hej = (ile/2) - 2 + spec_matrix[m][mm].move[k] + len(spectrum[mm].series)
                if (ile/2) - 2 + spec_matrix[m][mm].move[k] + len(spectrum[mm].series) > n and len(spec_matrix[m][mm].move) > 0:
                    del spec_matrix[m][mm].move[k:len(spec_matrix[m][mm].move)]
                    del spec_matrix[m][mm].rest[k:len(spec_matrix[m][mm].rest)]
                    # ile2 += 1
                    break

    # print("ile usu: " + str(ile2))
    # print("n: " + str(n))
    return spec_matrix


def complementary():
    global missing
    global spectrum2
    for o in range(len(spectrum2)):
        if spectrum2[o].komplementarny is False:
            for j in range(len(spectrum2)):
                if spectrum2[j].komplementarny is False and len(spectrum2[o].series) == len(spectrum2[j].series):
                    komp = True
                    for k in range(len(spectrum2[o].series)):
                        if spectrum2[o].series[k] == 'A':
                            if spectrum2[j].series[len(spectrum2[j].series) - 1 - k] != 'T':
                                komp = False
                        elif spectrum2[o].series[k] == 'T':
                            if spectrum2[j].series[len(spectrum2[j].series) - 1 - k] != 'A':
                                komp = False
                        elif spectrum2[o].series[k] == 'C':
                            if spectrum2[j].series[len(spectrum2[j].series) - 1 - k] != 'G':
                                komp = False
                        elif spectrum2[o].series[k] == 'G':
                            if spectrum2[j].series[len(spectrum2[j].series) - 1 - k] != 'C':
                                komp = False
                    if komp is True:
                        spectrum2[o].komplementarny = True
                        spectrum2[j].komplementarny = True
                        if spectrum2[o].min != spectrum2[j].min:
                            if spectrum2[o].min == -1 or spectrum2[j].min == -1:
                                spectrum2[o].min = -1
                                spectrum2[o].max = -1
                                spectrum2[o].max2 = -1
                                spectrum2[j].min = -1
                                spectrum2[j].max = -1
                                spectrum2[j].max2 = -1
                            elif spectrum2[o].min > spectrum2[j].min:
                                spectrum2[j].min = spectrum2[o].min
                                spectrum2[j].max = spectrum2[o].max
                                spectrum2[j].max2 = spectrum2[o].max2
                            else:
                                spectrum2[o].min = spectrum2[j].min
                                spectrum2[o].max = spectrum2[j].max
                                spectrum2[o].max2 = spectrum2[j].max2
                        break
    count = 0
    spectrum_size = len(spectrum2)
    for o in range(len(spectrum2)):
        if spectrum2[o].komplementarny is False:
            missing -= 1
            new_series = []
            for i in range(len(spectrum2[o].series)):
                new_series.append(complementary_sign(spectrum2[o].series[i]))
            #oli = Oligonucleotide()
            oli = copy.copy(spectrum2[o])
            new_series = new_series[::-1]
            oli.series = ''.join(new_series)
            oli.id = spectrum_size + count
            count += 1
            oli.komplementarny = True
            spectrum2[o].komplementarny = True
            spectrum2.append(copy.copy(oli))