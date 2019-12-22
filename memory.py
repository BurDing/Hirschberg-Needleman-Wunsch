import nose.tools as Test
import numpy as np
import argparse
from memory_profiler import memory_usage
import time

def get_delta():
    keys = ['A', 'C', 'T', 'G', '-']
    delta = {}
    for key in keys:
        delta[key] = {}
        for other in keys:
            if key == other:
                delta[key][other] = 1
            elif key == '-' or other == '-':
                delta[key][other] = -1
            else:
                delta[key][other] = -2
    delta['A']['G'] = -1
    delta['G']['A'] = -1
    delta['C']['T'] = -1
    delta['T']['C'] = -1
    delta['-']['-'] = float('-inf')

    return delta

def global_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
    computed by traceback_global.

    :param: v
    :param: w
    :param: delta
    """
    UP = (-1,0)
    LEFT = (0, -1)
    TOPLEFT = (-1, -1)
    ORIGIN = (0, 0)

    def traceback_global(v, w, pointers):
        i,j = len(v), len(w)
        new_v = []
        new_w = []
        while True:
            di, dj = pointers[i][j]
            if (di,dj) == LEFT:
                new_v.append('-')
                new_w.append(w[j-1])
            elif (di,dj) == UP:
                new_v.append(v[i-1])
                new_w.append('-')
            elif (di,dj) == TOPLEFT:
                new_v.append(v[i-1])
                new_w.append(w[j-1])
            i, j = i + di, j + dj
            if (i <= 0 and j <= 0):
                break
        return ''.join(new_v[::-1]), ''.join(new_w[::-1])

    M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]

    for i in range(len(v) + 1):
        for j in range(len(w) + 1):
            if i == 0 and j == 0:
                M[i][j] = 0
                pointers[i][j] = ORIGIN
            elif i == 0:
                M[i][j] = M[i][j-1] + delta['-'][w[j-1]]
                pointers[i][j] = LEFT
            elif j == 0:
                M[i][j] = M[i-1][j] + delta[v[i-1]]['-']
                pointers[i][j] = UP
            else:
                up_score = M[i-1][j] + delta[v[i-1]]['-']
                left_score = M[i][j-1] + delta['-'][w[j-1]]
                topleft_score = M[i-1][j-1] + delta[v[i-1]][w[j-1]]

                if up_score >= left_score and up_score >= topleft_score:
                    M[i][j] = up_score
                    pointers[i][j] = UP
                elif left_score > up_score and left_score >= topleft_score:
                    M[i][j] = left_score
                    pointers[i][j] = LEFT
                else:
                    M[i][j] = topleft_score
                    pointers[i][j] = TOPLEFT

    score = M[len(v)][len(w)]
    alignment_v, alignment_w = traceback_global(v, w, pointers)

    return score, alignment_v, alignment_w

def weight(a, b, delta):
        M = [[0 for i in range(2)] for j in range(len(a) + 1)]

        for j in range(len(b) + 1):
            for i in range(len(a) + 1):
                if i == 0 and j == 0:
                    M[0][0] = 0
                elif i == 0:
                    M[i][j%2] = M[i][(j-1)%2] + delta['-'][b[j-1]]
                elif j == 0:
                    M[i][j%2] = M[i-1][j%2] + delta[a[i-1]]['-']
                else:
                    up_score = M[i-1][j%2] + delta[a[i-1]]['-']
                    left_score = M[i][(j-1)%2] + delta['-'][b[j-1]]
                    topleft_score = M[i-1][(j-1)%2] + delta[a[i-1]][b[j-1]]

                    if up_score >= left_score and up_score >= topleft_score:
                        M[i][j%2] = up_score
                    elif left_score > up_score and left_score >= topleft_score:
                        M[i][j%2] = left_score
                    else:
                        M[i][j%2] = topleft_score

        return np.array(M)[:,len(b)%2]

def hirschberg(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment

    :param: v
    :param: w
    :param: delta
    """
    score, alignment_v, alignment_w = None, None, None

    if len(v) * len(w) == 0 or len(v) == 1 or len(w) == 1:
        score, alignment_v, alignment_w = global_align(v, w, delta)
    else:
        w_mid = len(w) // 2
        score_left = weight(v, w[:w_mid], delta)
        score_right = weight(v[::-1], w[w_mid:][::-1], delta)[::-1]

        v_mid = np.argmax(score_left + score_right)
        score_1, alignment_v_1, alignment_w_1 = hirschberg(v[:v_mid], w[:w_mid], delta)
        score_2, alignment_v_2, alignment_w_2 = hirschberg(v[v_mid:], w[w_mid:], delta)

        score = score_1 + score_2
        alignment_v = alignment_v_1 + alignment_v_2
        alignment_w = alignment_w_1 + alignment_w_2

    return score, alignment_v, alignment_w

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--length', type=int, default='100')
    parser.add_argument('--method', type=str, default='global')
    args = parser.parse_args()

    delta = get_delta()
    hiv1 = open('HIV1-Genome.txt', 'r').read()
    hiv2 = open('HIV2-Genome.txt', 'r').read()

    memo = 0
    if args.method == 'global':
        memo = max(memory_usage((global_align, (hiv1[0:args.length], hiv2[0:args.length], delta))))
    elif args.method == 'hirschberg':
        memo = max(memory_usage((hirschberg, (hiv1[0:args.length], hiv2[0:args.length], delta))))

    with open('result_' + args.method + '.txt', 'a') as f:
        f.write("%f," % memo)

if __name__ == '__main__':
    main()
