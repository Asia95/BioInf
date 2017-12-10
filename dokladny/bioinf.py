import string
import random
import math
import copy


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


class Seq:
    def __init__(self):
        self.chain = ""
        self.correct = False


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


spectrum2 = []
n = 0
sequence = ""
first_seq = ""
missing = 0
temp = 30


def spectrum():
    global spectrum2
    tmp = 30
    tmp2 = 32
    num_error = 1
    # num_repeated = 1
    plik = open('seq1.txt')
    try:
        f = plik.readlines()
        global first_seq
        first_seq = f[0]
        global sequence
        sequence = f[1]
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
                tmp_id = spectrum2[j].get_id() + 1
            if elem_spectrum[i] == first_seq.strip():
                f = 1
            else:
                f = 0
            o = Oligonucleotide(tmp_id, elem_spectrum[i], len(elem_spectrum[i]), 1, 0, f)
            spectrum2.append(o)

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
                choosen = random.radiant(0,len(spectrum2))
                # print(str(choosen))
                if spectrum2[choosen].min != 0:
                    # print(str(choosen))
                    spectrum2[choosen].min -= 1
                    ok = True

        i = 0
        # print(first_seq)
        while i < len(spectrum2):

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

            if spectrum2[i].min == 1 or spectrum2[i].min == 2:
                spectrum2[i].min = 1
                spectrum2[i].max = 2
                spectrum2[i].max2 = 4
            elif spectrum2[i].min == 3 or spectrum2[i].min == 4:
                spectrum2[i].min = 3
                spectrum2[i].max = 4
                spectrum2[i].max2 = -1
            elif spectrum2[i].min == 0:
                del spectrum2[i]
                i = i - 1
            else:
                spectrum2[i].min = -1
                spectrum2[i].max = -1
                spectrum2[i].max2 = -1

            i = i + 1

        return True
    else:
        return False

    for p in spectrum2: print(p)

bledy = ['A', 'T', 'C', 'G']

def find_first(spectrum):
    for o in range(len(spectrum)):
        if spectrum[o].first == 1:
            return spectrum[o]

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


def series_len(series):
    length = 0
    for s in series:
        if s != 'X':
            length += 1
    return length

# return index of complementary oligonucleotide or -1 if not found
def get_complementary(spectrum, id):
    tmp = ""
    for o in range(len(spectrum[id].series)):
        tmp += (complementary_sign(spectrum[id].series[o]))

    tmp = tmp[::-1]
    for o in range(len(spectrum)):
        if tmp == spectrum[o].series:
            return o
    return -1


def add_series(i, count, spectrum, series):
    tmp = spectrum[i].series[0:count]
    for t in range(len(tmp)):
        series.append(tmp[t])
    return series


def delete_series(i, count, spectrum, series):
    tmp = spectrum[i].series[0:count]
    for t in range(len(tmp)):
        if len(series) > 0:
            del series[-1]
    return series


def add_series_end(i, j, spectrum, series):
    tmp = spectrum[j].series
    for t in range(len(tmp)):
        series.insert(i, tmp[t])
        i += 1
    return series


def delete_series_end(i, j, spectrum, series):
    tmp = spectrum[j].series
    # for t in range(len(tmp)):
    del series[i:i + len(tmp)]
    # i += 1
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
    global spectrum2
    occur = []
    position = []
    new_seq = Seq()
    new_seq_list = []
    # sec_seq = []
    found = False
    count = 0
    # seq2 = ""
    global count2

    for o in range(len(spectrum)):
        if spectrum[o].times_used < spectrum[o].min:
            return False

    for o in range(len(seq)):
        if seq[o] == 'X':
            count += 1
            position.append(o)

    # for s in range(len(seq)):
        # seq2 += seq[s]

    # wp = 0
    # wz = 0
    # tmp = list(seq)
    new_seq = ""
    for i in range(len(seq)):
        # if seq[i] != 'X' and seq[i] != 'Y':
        if seq[i] != 'X':
            new_seq += seq[i]
    # new_seq = list(tmp)

    if len(new_seq) != n:
        new_seq_correct = False
    else:
        new_seq_correct = True

    join_chain = ''.join(new_seq)

    for s in range(len(spectrum2)):
        occur.clear()
        id1 = join_chain.find(spectrum2[s].series)
        # id1 = new_seq_list[i].chain.index(spectrum2[s].series)
        if id1 != -1:
            occur.append(id1)
        while id1 >= 0:
            # join_chain = ''.join(new_seq)
            id1 = join_chain.find(spectrum2[s].series, id1 + 1)
            # id1 = new_seq_list[i].chain.find(spectrum2[s].series, id1+1)
            if id1 != -1:
                occur.append(id1)

        if spectrum2[s].min == 0:
            if len(occur) > 0 and len(occur) < 2:
                new_seq_correct = False
        elif spectrum2[s].min == 1:
            if len(occur) > 1 and len(occur) < 4:
                new_seq_correct = False
        elif spectrum2[s].min == 3:
            if len(occur) > 3:
                new_seq_correct = False
        elif spectrum2[s].min == -1:
            if len(occur) >= 4:
                new_seq_correct = False

    if new_seq_correct == True:
        count2 += 1
        print(str(count2))
        print("wygenerowana sekwencja: " + join_chain)
        print("jakość dopasowania: " + '%.2f' % quality(seq) + " %")
        print()
        found = True
    # count2 += 1
    # if count2 == 1:
    # print("hej")
    if found == True:
        return True
    else:
        return False


def get_id(id2, spectrum2):
    for o in range(len(spectrum2)):
        if spectrum2[o].id == id2:
            return o
    return -1


def odciecie(spec_matrix, spectrum, seq, accual_id):
    min_lenght = 1000
    # missing3 = 0
    missing1 = 0
    missing2 = 0

    for s in range(len(spectrum)):
        if spectrum[s].times_used < spectrum[s].min:
            min_shift = n
            for o in range(len(spectrum)):
                if spectrum[o].times_used < spectrum[o].max2:
                    if len(spec_matrix[s][o].move) > 0:
                        if s == o and len(spec_matrix[s][o].move) > 1 and spec_matrix[s][o].move[1] < min_shift:
                            min_shift = spec_matrix[s][o].move[1]
                        elif spec_matrix[s][o].move[0] == 0 and len(spec_matrix[s][o].move) > 1 and spec_matrix[s][o].move[1] < min_shift:
                            min_shift = spec_matrix[s][o].move[1]
                        elif spec_matrix[s][o].move[0] < min_shift:
                            min_shift = spec_matrix[s][o].move[0]
                    else:
                        min_shift = 1

                    if min_shift == 1:
                        break

            if temp_series(spectrum[s].series) == temp:
                if min_shift < n:
                    missing1 += (spectrum[s].min - spectrum[s].times_used) * min_shift
                    if len(spectrum[s].series) < min_lenght:
                        min_lenght = len(spectrum[s].series)
                else:
                    missing1 += spectrum[s].min - spectrum[s].times_used
                    if len(spectrum[s].series) < min_lenght:
                        min_lenght = len(spectrum[s].series)
            else:
                if min_shift < n:
                    missing2 += (spectrum[s].min - spectrum[s].times_used) * min_shift
                    if len(spectrum[s].series) < min_lenght:
                        min_lenght = len(spectrum[s].series)
                else:
                    missing2 += spectrum[s].min - spectrum[s].times_used
                    if len(spectrum[s].series) < min_lenght:
                        min_lenght = len(spectrum[s].series)

    if missing2 > missing1:
        missing3 = missing2
    else:
        missing3 = missing1

    if missing3 > 0:
        if len(spectrum[accual_id].series) < (missing3 / 2) - 1 + min_lenght - 1:
            if series_len(seq) + 1 + missing3 / 2 - 1 + min_lenght - 1 > n:
                return True
            else:
                return False
        else:
            if series_len(seq) + len(spectrum[accual_id].series) > n:
                return True
            else:
                return False
    else:
        return False

def odciecie2(spec_matrix, spectrum, seq, accual_id):
    min_lenght = 10000
    missing1 = 0
    missing2 = 0

    for s in range(len(spectrum)):
        if spectrum[s].times_used < spectrum[s].min:
            if temp_series(spectrum[s].series) == temp:
                missing1 += (spectrum[s].min - spectrum[s].times_used)
                if len(spectrum[s].series) < min_lenght:
                    min_lenght = len(spectrum[s].series)
            else:
                missing2 += (spectrum[s].min - spectrum[s].times_used)
                if len(spectrum[s].series) < min_lenght:
                    min_lenght = len(spectrum[s].series)

    if missing2 > missing1:
        missing3 = missing2
    else:
        missing3 = missing1

    if missing3 > 0:
        if len(spectrum[accual_id].series) < missing3/2 - 1 + min_lenght - 1:
            if series_len(seq) + 1 + missing3/2 - 1 + min_lenght - 1 > n:
                return True
            else:
                return False
        else:
            if series_len(seq) + len(spectrum[accual_id].series) > n:
                return True
            else:
                return False
    else:
        return False

def set_available(id1, id2, tmp_spectrum, tmp_spec_matrix):
    if tmp_spectrum[id1].max2 != -1:
        if tmp_spectrum[id1].times_used >= tmp_spectrum[id1].max2:
            for s in range(len(tmp_spec_matrix)):
                tmp_spec_matrix[s][id1].available = False
                tmp_spec_matrix[s][id2].available = False
        else:
            for s in range(len(tmp_spec_matrix)):
                tmp_spec_matrix[s][id1].available = True
                tmp_spec_matrix[s][id2].available = True
    return tmp_spec_matrix


def del_times_used(spectrum, i):
    if (i == 1 or i == 0) and spectrum[i].times_used == 1:
        print()
    else:
        spectrum[i].times_used -= 1
    return spectrum


def find_next(actual, spec_matrix, spectrum, seq, it, found, missing):
    i = it
    ifcopy = False

    if spec_matrix[actual][i].available == True:
        if (spectrum[actual].id != spectrum[i].id) or (spectrum[actual].id == spectrum[i].id and spectrum[actual].max2 == -1) or (spectrum[actual].id == spectrum[i].id and spectrum[actual].times_used < spectrum[actual].max2 - 1):
            for j in range(len(spec_matrix[actual][i].move)):
                if spec_matrix[actual][i].rest[j] == 0:
                    if series_len(seq) + spec_matrix[actual][i].move[j] <= n:
                        if ifcopy == True:
                            graph(actual, i, j, tmp_spec_matrix, tmp_spectrum, seq, missing)
                        else:
                            graph(actual, i, j, spec_matrix, spectrum, seq, missing)
                        found = True
                else:
                    ifcopy = True
                    tmp_spectrum = copy.deepcopy(spectrum)
                    tmp_spec_matrix = copy.deepcopy(spec_matrix)
                    tmp_spectrum[i].times_used += 1
                    i2 = get_complementary(tmp_spectrum, i)
                    tmp_spectrum[i2].times_used += 1
                    tmp_spec_matrix = set_available(i2, i, tmp_spectrum, tmp_spec_matrix)
    return found


# count3 = 0
def graph(before, actual, it, spec_matrix, spectrum, seq, missing):

    # print("actual : " + str(actual))
    # print("before : " + str(before))

    # plik = open('plik.txt', 'a')
    # plik.write(''.join(seq) + "\n")
    # plik.write(str(before) + "\n")
    # plik.close()

    # global count3

    # if actual == 55 and before == 2:
        # count3 += 1
        # if count3 == 1:
            # print("stop")

    # count3 += 1
    # if count3 == 3793:
        # print("stop")

    found = False
    not_ok = False
    added_more = False
    added_one = False
    added_x = False
    index1 = -1
    index2 = -1
    if before != -1:
        seq_len1 = len(seq)
        seq = add_series(before, spec_matrix[before][actual].move[it], spectrum, seq)
        # if len(seq) > 8:
            # if seq[8] == 'T':
                # print("16")
        seq_len2 = len(seq)
        if seq_len2 - seq_len1 > 1:
            added_more = True
        else:
            if seq_len2 - seq_len1 == 1:
                added_one = True
        missing += (spec_matrix[before][actual].move[it] - 1) * 2
        if spec_matrix[before][actual].move[it] == len(spectrum[before].series):
            seq.append('X')
            added_x = True
        index1 = get_complementary(spectrum, actual)
        spectrum[actual].times_used += 1
        spectrum[index1].times_used += 1
        spec_matrix = set_available(index1, actual, spectrum, spec_matrix)

    if series_len(seq) + len(spectrum[actual].series) <= n:
        if odciecie2(spec_matrix, spectrum, seq, actual) is False:
            for i in range(len(spec_matrix[actual])):
                found = find_next(actual, spec_matrix, spectrum, seq, i, found, missing)
    else:
        frag = spectrum[before].series[spec_matrix[before][actual].move[it]:spec_matrix[before][actual].move[it] + len(spectrum[before].series)]
        seq += frag
        if check_if_ok(spectrum, seq) is False:
            not_ok = True
            # print("Wygenerował niepoprawną sekwencję1: ")
            # for i in range(len(seq)):
            # print(seq)
        else:
            not_ok = True
    if series_len(seq) + len(spectrum[actual].series) == n:
        tmp2 = len(seq)
        seq = add_series_end(tmp2, actual, spectrum, seq)
        index2 = get_complementary(spectrum, actual)
        spectrum[actual].times_used += 1
        spectrum[index2].times_used += 1

        if check_if_ok(spectrum, seq) is False:
            # print("Wygenerował niepoprawną sekwencję2: ")
            # print(seq)
            seq = delete_series_end(tmp2, actual, spectrum, seq)
            if added_x == True:
                added_x = False
                if len(seq) > 0:
                    del seq[-1]
            if index2 != -1:
                spectrum = del_times_used(spectrum, index2)
                spectrum = del_times_used(spectrum, actual)
        else:
            seq = delete_series_end(tmp2, actual, spectrum, seq)
            if added_x is True:
                added_x = False
                if len(seq) > 0:
                    del seq[-1]
            if index2 != -1:
                spectrum = del_times_used(spectrum, index2)
                spectrum = del_times_used(spectrum, actual)

    elif not_ok is True:
        not_ok = False
        for i in range(len(frag)):
            if len(seq) > 0:
                del seq[-1]
        if added_x is True:
            added_x = False
            if len(seq) > 0:
                del seq[-1]
    if added_one is True:
        added_one = False
        if len(seq) > 0:
            del seq[-1]
        if added_x is True:
            added_x = False
            if len(seq) > 0:
                del seq[-1]
    if added_more is True:
        added_more = False
        seq = delete_series(before, spec_matrix[before][actual].move[it], spectrum, seq)
        if added_x is True:
            added_x = False
            if len(seq) > 0:
                del seq[-1]
    if index1 != -1:
        spectrum = del_times_used(spectrum, index1)
        spectrum = del_times_used(spectrum, actual)
        spec_matrix = set_available(index1, actual, spectrum, spec_matrix)


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
            new_series = ""
            for i in range(len(spectrum2)):
                new_series.append(complementary_sign(spectrum2[o].series[i]))
            # oli = Oligonucleotide()
            oli = spectrum2[o]
            new_series = new_series[::-1]
            oli.series = new_series
            oli.id = spectrum_size + count
            count += 1
            oli.komplementarny = True
            spectrum2[o].komplementarny = True
            spectrum2.append(list(oli))


if __name__ == '__main__':

    spectrum()
    complementary()
    # spectrum3 = spectrum2
    ma = matrix(spectrum2)

    oli_first = find_first(spectrum2)
    actual = get_id(oli_first.id, spectrum2)
    i = get_complementary(spectrum2, actual)
    spectrum2[actual].times_used += 1
    spectrum2[i].times_used += 1
    ma = set_available(i, actual, spectrum2, ma)
    seq = []
    graph(-1, actual, -1, ma, spectrum2, seq, 0)
    print("koniec")
