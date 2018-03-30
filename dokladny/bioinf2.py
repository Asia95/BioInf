import copy
import time
import math
import random
from spectrum2 import *


class Seq:
    def __init__(self):
        self.chain = ""
        self.correct = False


missing = 0

bledy = ['A', 'T', 'C', 'G']


def find_first(graph):
    for o in range(len(graph)):
        if graph[o].start.first == 1:
            return graph[o].start


def series_len(series):
    length = 0
    for s in series:
        if s != 'X':
            length += 1
    return length

# return index of complementary oligonucleotide or -1 if not found
def get_complementary(spectrum, id):
    tmp = ""
    for o in range(len(spectrum[id].start.series)):
        tmp += (complementary_sign(spectrum[id].start.series[o]))

    tmp = tmp[::-1]
    for o in range(len(spectrum)):
        if tmp == spectrum[o].start.series:
            return o
    return -1


def add_series(i, count, spectrum, series):
    tmp = spectrum[i].start.series[0:count]
    for t in range(len(tmp)):
        series.append(tmp[t])
    return series


def delete_series(i, count, spectrum, series):
    tmp = spectrum[i].start.series[0:count]
    for t in range(len(tmp)):
        if len(series) > 0:
            del series[-1]
    return series


def add_series_end(i, j, spectrum, series):
    tmp = spectrum[j].start.series
    for t in range(len(tmp)):
        series.insert(i, tmp[t])
        i += 1
    return series


def delete_series_end(i, j, spectrum, series):
    tmp = spectrum[j].start.series
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
    occur = []
    position = []
    found = False
    count = 0
    global count2

    for o in range(len(spectrum)):
        if spectrum[o].start.times_used < spectrum[o].start.min and spectrum[o].start.times_used > spectrum[o].start.max:
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
        if seq[i] != 'X':
            new_seq += seq[i]
    # new_seq = list(tmp)

    if len(new_seq) != n:
        new_seq_correct = False
    else:
        new_seq_correct = True

    join_chain = ''.join(new_seq)

    for s in range(len(graph)):
        occur.clear()
        id1 = join_chain.find(graph[s].start.series)
        # id1 = new_seq_list[i].chain.index(spectrum2[s].series)
        if id1 != -1:
            occur.append(id1)
        while id1 >= 0:
            # join_chain = ''.join(new_seq)
            id1 = join_chain.find(graph[s].start.series, id1 + 1)
            # id1 = new_seq_list[i].chain.find(spectrum2[s].series, id1+1)
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
        print("Generated sequence: " + join_chain)
        print("Quality of the match: " + '%.2f' % quality(seq) + " %")
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
        if spectrum2[o].start.id == id2:
            return o
    return -1


def odciecie2(graph, seq, accual_id):
    min_lenght = 10000
    missing1 = 0
    missing2 = 0

    for s in range(len(graph)):
        if graph[s].start.times_used < graph[s].start.min:
            if temp_series(graph[s].start.series) == temp:
                missing1 += (graph[s].start.min - graph[s].start.times_used)
                if len(graph[s].start.series) < min_lenght:
                    min_lenght = len(graph[s].start.series)
            else:
                missing2 += (graph[s].start.min - graph[s].start.times_used)
                if len(graph[s].start.series) < min_lenght:
                    min_lenght = len(graph[s].start.series)

    if missing2 > missing1:
        missing3 = missing2
    else:
        missing3 = missing1

    if missing3 > 0:
        if len(graph[accual_id].start.series) < missing3/2 - 1 + min_lenght - 1:
            if series_len(seq) + 1 + missing3/2 - 1 + min_lenght - 1 > n:
                return True
            else:
                return False
        else:
            if series_len(seq) + len(graph[accual_id].start.series) > n:
                return True
            else:
                return False
    else:
        return False

def set_available(id1, id2, graph):
    if graph[id1].start.max2 != -1:
        if graph[id1].start.times_used >= graph[id1].start.max2:
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


def del_times_used(spectrum, i):
    if (i == 1 or i == 0) and spectrum[i].start.times_used == 1:
        print()
    else:
        spectrum[i].start.times_used -= 1
    return spectrum


def find_next(actual, graph, seq, it, found, missing):
    i = it
    ifcopy = False

    if graph[actual].edges[i].available == True:
        idx = get_id(graph[actual].edges[i].end.id, graph)
        if (graph[actual].start.id != graph[idx].start.id) or (graph[actual].start.id == graph[idx].start.id and graph[actual].start.max2 == -1) or (graph[actual].start.id == graph[idx].start.id and graph[actual].start.times_used < graph[actual].start.max2 - 1):
            for j in range(len(graph[actual].edges[i].move)):
                if graph[actual].edges[i].rest[j] == 0:
                    if series_len(seq) + graph[actual].edges[i].move[j] <= n:
                        if ifcopy == True:
                            find_seq(actual, idx, j, tmp_graph, seq, missing)
                        else:
                            find_seq(actual, idx, j, graph, seq, missing)
                        found = True
                else:
                    ifcopy = True
                    tmp_graph = copy.deepcopy(graph)
                    tmp_graph[idx].start.times_used += 1
                    i2 = get_complementary(tmp_graph, idx)
                    tmp_graph[i2].start.times_used += 1
                    tmp_graph = set_available(i2, idx, tmp_graph)
    return found


# count3 = 0
def find_seq(before, actual, it, graph, seq, missing):

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

    for i in range(len(graph[before].edges)):
        if actual == graph[before].edges[i].end.id:
            actual_idx = i
            break


    found = False
    not_ok = False
    added_more = False
    added_one = False
    added_x = False
    index1 = -1
    index2 = -1
    if before != -1:

        seq_len1 = len(seq)
        seq = add_series(before, graph[before].edges[actual_idx].move[it], graph, seq)
        # if len(seq) > 8:
            # if seq[8] == 'T':
                # print("16")
        seq_len2 = len(seq)
        if seq_len2 - seq_len1 > 1:
            added_more = True
        else:
            if seq_len2 - seq_len1 == 1:
                added_one = True
        missing += (graph[before].edges[actual_idx].move[it] - 1) * 2
        if graph[before].edges[actual_idx].move[it] == len(graph[before].start.series):
            seq.append('X')
            added_x = True
        index1 = get_complementary(graph, actual)
        graph[actual].start.times_used += 1
        graph[index1].start.times_used += 1
        graph = set_available(index1, actual, graph)

    if series_len(seq) + len(graph[actual].start.series) <= n:
        if odciecie2(graph, seq, actual) is False:
            for i in range(len(graph[actual].edges)):
                found = find_next(actual, graph, seq, i, found, missing)
    else:
        frag = graph[before].start.series[graph[before].edges[actual_idx].move[it]:graph[before].edges[actual_idx].move[it] + len(graph[before].start.series)]
        seq += frag
        if check_if_ok(graph, seq) is False:
            not_ok = True
            # print("Wygenerował niepoprawną sekwencję1: ")
            # for i in range(len(seq)):
            # print(seq)
        else:
            not_ok = True
    if series_len(seq) + len(graph[actual].start.series) == n:
        tmp2 = len(seq)
        seq = add_series_end(tmp2, actual, graph, seq)
        index2 = get_complementary(graph, actual)
        graph[actual].start.times_used += 1
        graph[index2].start.times_used += 1

        if check_if_ok(graph, seq) is False:
            # print("Wygenerował niepoprawną sekwencję2: ")
            # print(seq)
            seq = delete_series_end(tmp2, actual, graph, seq)
            if added_x == True:
                added_x = False
                if len(seq) > 0:
                    del seq[-1]
            if index2 != -1:
                graph = del_times_used(graph, index2)
                graph = del_times_used(graph, actual)
        else:
            seq = delete_series_end(tmp2, actual, graph, seq)
            if added_x is True:
                added_x = False
                if len(seq) > 0:
                    del seq[-1]
            if index2 != -1:
                graph = del_times_used(graph, index2)
                graph = del_times_used(graph, actual)

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
        seq = delete_series(before, graph[before].edges[actual_idx].move[it], graph, seq)
        if added_x is True:
            added_x = False
            if len(seq) > 0:
                del seq[-1]
    if index1 != -1:
        graph = del_times_used(graph, index1)
        graph = del_times_used(graph, actual)
        graph = set_available(index1, actual, graph)


if __name__ == '__main__':

    spectrum()
    n = return_len_seq()
    sequence = return_seq()
    complementary()
    matrix()

    start = time.clock()

    oli_first = find_first(graph)
    actual = get_id(oli_first.id, graph)
    i = get_complementary(graph, actual)
    graph[actual].start.times_used += 1
    graph[i].start.times_used += 1
    graph = set_available(i, actual, graph)
    seq = []
    find_seq(-1, actual, -1, graph, seq, 0)

    print ("time: " + str(time.clock() - start))

    print("end")
