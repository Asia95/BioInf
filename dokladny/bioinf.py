import copy
import time
from spectrum import *


class Seq:
    def __init__(self):
        self.chain = ""
        self.correct = False


missing = 0

bledy = ['A', 'T', 'C', 'G']

def find_first(spectrum):
    for o in range(len(spectrum)):
        if spectrum[o].first == 1:
            return spectrum[o]


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


misteake = ['A','T','C','G']
per = []
per2 = []


def make_per2(count, position, comb):
    tmp = ""
    for i in range(len(misteake)):
        comb[position] = misteake[i]
        if position == count - 1:
            for j in range(len(comb)):
                tmp += comb[j]
            per2.append(tmp)
            tmp = ""
        else:
            make_per2(count, position+1, comb)


def make_per(count, used, countSpace, countSign, road):
    if count == countSpace - 1:
        road[count] = countSign - used
        per.append(road)
    else:
        for i in range(countSign - used):
            road[count] = i
            make_per(count+1, used+1, countSpace, countSign, road)


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
        if spectrum[o].times_used < spectrum[o].min or spectrum[o].times_used > spectrum[o].max2:
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


    if count != 0:
        road = []
        for i in range(count):
            road.append(0)
        comb = []
        for i in range(n-(len(seq)-count)):
            comb.append(0)
        if n-(len(seq)-count-(n-(len(seq)-count))) > 0:
            make_per(0, 0, count, n - len(seq) - count - (n - len(seq) - count), road)
            make_per2(n - len(seq) - count - (n - len(seq) - count), 0, comb)

    seq2 = ""
    for i in (range(len(seq))):
        seq2 += seq[i]

    seqtmp = seq2
    for i in range(len(per)):
        for j in range(len(per2)):
            seqtmp = seq2
            lenPos = len(position) - 1
            lenPer2 = len(per2[j])
            k = len(seqtmp)
            while k >= 0:
                if k == position[lenPos]:
                    seqtmp.pop()
                    seqtmp.insert(k, per2[j][lenPer2 - per[i][lenPos] : lenPer2 - per[i][lenPos] + per[i][lenPos]])
                    if lenPos > 0:
                        lenPer2 = lenPer2 - per[i][lenPos]
                        lenPos -= 1
                    else:
                        new_seq.chain = seqtmp
                        new_seq.correct = True
                        new_seq_list.append(new_seq)
                        break
                k -= 1

    if len(new_seq_list) == 0:
        seqtmp = ""
        for i in range(len(seq2)):
            if seq2[i] != 'X':
                seqtmp += seq2[i]
        seq2 = seqtmp
        new_seq.chain = seq2
        new_seq.correct = True
        new_seq_list.append(new_seq)
        # new_seq = list(tmp)

    for i in range(len(new_seq_list)):
        if len(new_seq_list[i].chain) != n:
            new_seq_list[i].correct = False
        else:
            new_seq_list[i].correct = True

        join_chain = ''.join(new_seq_list[i].chain)

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
                if len(occur) > 0 and len(occur) < 1:
                    new_seq_list[i].correct = False
            elif spectrum2[s].min == 1:
                if len(occur) > 1 and len(occur) < 3:
                    new_seq_list[i].correct = False
            elif spectrum2[s].min == 2:
                if len(occur) > 2 and len(occur) < 5:
                    new_seq_list[i].correct = False
            elif spectrum2[s].min == 4:
                if len(occur) > 4:
                    new_seq_list[i].correct = False
            elif spectrum2[s].min == -1:
                if len(occur) >= 5:
                    new_seq_list[i].correct = False

        if new_seq_list[i].correct == True:
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


def check_correct(spec_matrix, spectrum, seq, accual_id):
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
        if check_correct(spec_matrix, spectrum, seq, actual) is False:
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


if __name__ == '__main__':

    spectrum()
    n = return_len_seq()
    sequence = return_seq()
    complementary()
    ma = matrix(spectrum2)

    start = time.clock()

    oli_first = find_first(spectrum2)
    actual = get_id(oli_first.id, spectrum2)
    i = get_complementary(spectrum2, actual)
    spectrum2[actual].times_used += 1
    spectrum2[i].times_used += 1
    ma = set_available(i, actual, spectrum2, ma)
    seq = []
    graph(-1, actual, -1, ma, spectrum2, seq, 0)

    print ("time: " + str(time.clock() - start))

    print("koniec")
