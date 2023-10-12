from Bio import SeqIO
import numpy as np
import subprocess
from tqdm import tqdm

def read_sam(path_to_reads):
    sequence = []
    with open(path_to_reads, "r") as sam_file:
        for line in sam_file:
            if line.startswith("@"):  # пропускаем заголовочные строки
                continue
            parts = line.strip().split("\t")
            sequence.append(parts[9])  # получаем последовательность
    return sequence

def read_fasta(path_to_reads):
    sequence = []
    with open(path_to_reads, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                sequence.append(line)
    return sequence

def read_fastq(path_to_reads):
    sequences = []
    with open(path_to_reads, "r") as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            sequence = lines[i+1].strip()
            sequences.append(sequence)
    return sequences



def write_fastq(file_path, sequence, quality):
    with open(file_path, "w") as file:
        file.write("@consensus\n")  # Заголовок последовательности
        file.write(sequence + "\n")     # Сама последовательность
        file.write("+\n")                # Заголовок качества
        file.write(quality + "\n")       # Качество для последовательности

def fastq_to_fasta(path_fastq, path_fasta):
    # Функция перевода fastq в fasta формат
    # На вход путь до fastq файла и название fasta файла

    with open(path_fasta, "w") as output:
        sequences = SeqIO.parse(path_fastq, "fastq")
        SeqIO.write(sequences, output, "fasta")
    return

def read_seq_from_file(path):
    # На вход принимается путь до fasta файла
    # На выходе получаем список с ридами
    sequences = []
    with open(path, 'r') as file:
        seq = ''
        for line in file:
            if line.startswith('>'):
                if seq != '':
                    sequences.append(seq)
                seq = ''
            else:
                seq += line.strip()
        sequences.append(seq)
    return sequences

def write_seq_in_file_with_length(name_file, list_seq, start, end):
    # Функция записи среза последовательности в файл
    # На вход подаётся название файла, список последовательностей, начало среза, конец среза
    # На выходе получаем файл со срезами последовательностей от start до end, end не включительно
    # name_file - string - название файла
    # list_seq - list - список последовательностей
    # start - integer - начало среза
    # end - integer - конец среза
    with open(name_file, "w") as f:
        for j, i in enumerate(list_seq):
            header = f">sequence_id{j}\n"
            f.write(header)
            if start == end:
                f.write(i +'\n')
            else:
                f.write(i[start:end] +'\n')
    return

def write_seq_in_file_with_length_and_name(name_file, list_seq, name_seq, start, end):
    # Функция записи среза последовательности в файл
    # На вход подаётся название файла, список последовательностей, начало среза, конец среза
    # На выходе получаем файл со срезами последовательностей от start до end, end не включительно
    # name_file - string - название файла
    # list_seq - list - список последовательностей
    # start - integer - начало среза
    # end - integer - конец среза
    with open(name_file, "w") as f:
        for j, i in enumerate(list_seq):
            header = f">sequence_id{name_seq[j]}\n"
            f.write(header)
            if start == end:
                f.write(i +'\n')
            else:
                f.write(i[start:end] +'\n')
    return

def muscle(path_in, path_out):
    # Функция запуска программы MUSCLE с штрафом за пропуск = 2
    command = f"./muscle3.8.31_i86linux64 -gapopen -2 -in {path_in} -out {path_out} -quiet"
    subprocess.run(command, shell=True)
    return


def muscle_with_gap1000(path_in, path_out):
    # Функция запуска программы MUSCLE с штрафом за пропуск = 1000
    command = f"./muscle3.8.31_i86linux64 -gapopen -2000 -in {path_in} -out {path_out} -quiet"
    subprocess.run(command, shell=True)
    return

#Процент каждого нуклеотида в столбце
def per_nucl_in_cal(allig):
    count = []
    for i in range(len(allig[0])):
        tmp = np.zeros(4)
        for j in range(len(allig)):
            for k,l in enumerate(['A','T','G','C']):
                if allig[j][i] == l:
                    tmp[k] += 1 / len(allig)
                    break
        count.append(tmp)
    return count

#Подсчета процента гэпов в столбце
def count_gaps (allig):
    count_gap = []
    for i in range (len(allig[0])):
        tmp = 0
        for j in allig:
            if j[i] == '-':
                tmp += 1/len(allig)
        count_gap.append(tmp)
    return (count_gap)


#Удалени столбцов, где больше 50% гэпов

def delete_gap (allig, count, count_nucleotide, perc_gap):
    for k,l in enumerate (allig):
        a = l
        for i, j in enumerate (count):
            if j >= perc_gap and count_nucleotide[i][0] <= 0.5 and count_nucleotide[i][1] <= 0.5 and count_nucleotide[i][2] <= 0.5 and count_nucleotide[i][3] <= 0.5:
                a = a[:i] + 'f' + a[i+1:]       
        allig[k] = a

    for i,j in enumerate(allig):
        allig[i] =j.replace('f','')
    return (allig)

#Вероятность находиться в состоянии в начальный момент времени
def start_probability (allig):
    tmp = 0
    pi = {}
    for i in allig:
        if i[0] != '-':
            tmp += 1
    if tmp >= len(allig)/2:
        pi['M'] = (tmp + 1)/(len(allig) + 3)
        pi['D'] = (len(allig) - tmp + 1)/(len(allig) + 3)
        pi['I'] = 1/(len(allig) + 3)
    else: 
        pi['M'] = 1/(len(allig) + 3)
        pi['D'] = 1/(len(allig) + 3)
        pi['I'] = len(allig)/(len(allig) + 3)
    return pi



#Буквенные переходы 3 - M, 2 - D, 1 - I
import numpy as np

def matrix_of_MDI (allig):
    a = np.ones((len(allig), len(allig[0]) - 1))

    for i in range(len(allig[0]) - 1):
        count1 = 0
        count2 = 0
        for j in range(len(allig)):
            if allig[j][i] == '-':
                count1 += 1
            if allig[j][i + 1] == '-':
                count2 += 1
        if (count1 >= len(allig)/2) and (count2 >= len(allig)/2):
            for j in range(len(allig)):
                if allig[j][i] == '-' and allig[j][i + 1] == '-' and i != 0:
                    a[j,i] = 11
                elif allig[j][i] == '-' and allig[j][i + 1] != '-':
                    a[j,i] = 13
                elif allig[j][i] != '-' and allig[j][i + 1] == '-':
                    a[j,i] = 11
                else:
                    a[j,i] = 11
        elif count2 >= len(allig)/2:
            for j in range(len(allig)):
                if allig[j][i] == '-' and allig[j][i + 1] == '-':
                    a[j,i] = 21
                elif allig[j][i] == '-' and allig[j][i + 1] != '-':
                    a[j,i] = 21
                elif allig[j][i] != '-' and allig[j][i + 1] == '-':
                    a[j,i] = 31
                else:
                    a[j,i] = 31
        elif count1 >= len(allig)/2:
            for j in range(len(allig)):
                if allig[j][i] == '-' and allig[j][i + 1] == '-' and i != 0:
                    a[j,i] = 12
                elif allig[j][i] == '-' and allig[j][i + 1] != '-':
                    a[j,i] = 13
                elif allig[j][i] != '-' and allig[j][i + 1] == '-':
                    a[j,i] = 12
                else:
                    a[j,i] = 13
        else:
            for j in range(len(allig)):
                if allig[j][i] == '-' and allig[j][i + 1] == '-':
                    a[j,i] = 22
                elif allig[j][i] == '-' and allig[j][i + 1] != '-':
                    a[j,i] = 23
                elif allig[j][i] != '-' and allig[j][i + 1] == '-':
                    a[j,i] = 32
                else:
                    a[j,i] = 33
    return a


# Вероятности переходов
def probability_of_trans(allig, matrix):
    A = {}
    for i in range(len(allig[0]) - 1):
        count1 = 0
        count2 = 0
        count3 = 0
        count11 = 0
        count22 = 0
        count33 = 0
        count13 = 0
        count23 = 0
        count31 = 0
        count32 = 0
        count21 = 0
        count12 = 0
        for j in range(len(allig)):
            if matrix[j][i] == 33:
                count3 += 1
                count33 += 1
            elif matrix[j][i] == 22:
                count2 += 1
                count22 += 1
            elif matrix[j][i] == 11:
                count1 += 1
                count11 += 1
            elif matrix[j][i] == 13:
                count1 += 1
                count13 += 1
            elif matrix[j][i] == 23:
                count2 += 1
                count23 += 1
            elif matrix[j][i] == 31:
                count3 += 1
                count31 += 1
            elif matrix[j][i] == 21:
                count2 += 1
                count21 += 1
            elif matrix[j][i] == 12:
                count1 += 1
                count12 += 1
            else:
                count3 += 1
                count32 += 1
        if count3 != 0:
            A[f'{i}MM'] = (count33 + 1)/(count3 + 3)
            A[f'{i}MI'] = (count31 + 1)/(count3 + 3)
            A[f'{i}MD'] = (count32 + 1)/(count3 + 3)
        if count2 != 0:
            A[f'{i}DI'] = (count21 + 1)/(count2 + 3)
            A[f'{i}DD'] = (count22 + 1)/(count2 + 3)
            A[f'{i}DM'] = (count23 + 1)/(count2 + 3)
        if count1 != 0:
            A[f'{i}ID'] = (count12 + 1)/(count1 + 3)
            A[f'{i}II'] = (count11 + 1)/(count1 + 3)
            A[f'{i}IM'] = (count13 + 1)/(count1 + 3)
    return A


#Выходные вероятности
sost = 'MDI'

def count_nucl_in_cal (allig):
    count = []
    for i in range (len(allig[0])):
        tmp = np.zeros(4)
        for j in range (len(allig)):
            for k,l in enumerate(['A','T','G','C']):
                if allig[j][i] == l:
                    tmp[k] += 1
                    break
        count.append(tmp)
    return count



def percetn_of_nucl_in_cal(count):
    nucl = ['A','T','G','C']
    per = {}
    for s in 'MDI':
        for i,k in enumerate (count):
            tmp = []
            for j,l in enumerate (k):
                if s == 'M':
                    per[f'{s}{i}{nucl[j]}'] = (l+1)/(sum(k) + 4)
                elif s == 'I':
                    per[f'{s}{i}{nucl[j]}'] = 0.25
                else:
                    per[f'{s}{i}{nucl[j]}'] = 1
    return per



def consensus(allig, out_pr):
    nucl = ['A','T','G','C']
    consensus = ''
    for i in range (len(allig[0])):
        tmp = []
        for j in nucl:
            tmp.append(out_pr[f'M{i}{j}'])
        consensus += nucl[tmp.index(max(tmp))] 
    return consensus



def delete_insert(alligment):
    list_nucl = []
    max_length = max(alligment, key=len)
    for i in range(len(max_length)):
        tmp_list = []
        for j in range(len(alligment)):
            try:
                tmp_list.append(alligment[j][i])
            except:
                continue
        list_nucl.append(tmp_list)
    list_seq_without_insert = []

    for i in list_nucl:
        if i.count('-') < int(len(i) / 2):
            list_seq_without_insert.append(i)
    list_seq_without_insert_row = []
    max_length = max(list_seq_without_insert, key=len)


    for i in range(len(max_length)):
        tmp_list = []
        for j in range(len(list_seq_without_insert)):
            try:
                tmp_list.append(list_seq_without_insert[j][i])
            except:
                continue
        list_seq_without_insert_row.append(tmp_list)
    for i in range(len(list_seq_without_insert_row)):
        for j in range(len(list_seq_without_insert_row[i]) - 1, -1, -1):
            if list_seq_without_insert_row[i][j] == '-':
                del list_seq_without_insert_row[i][j]
    okon_list = []

    for i in list_seq_without_insert_row:
        okon_list.append(''.join(i))
    return okon_list


def sort_fasta(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        records.append(record)
    sorted_records = sorted(records, key=lambda x: int(x.description.split("id")[1]))
    return sorted_records


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse = dna[::-1]
    reverse_complement = ''.join([complement[base] for base in reverse])
    return reverse_complement



def quality(list_seq):
    list_quality = []
    for i in range(len(list_seq[0][:-100])):
        count = 0
        for j in range(len(list_seq) - 1):
            if list_seq[0][i] == list_seq[j + 1][i]:
                count += 1
        list_quality.append(count / (len(list_seq) - 1))
    return list_quality

def convert_to_phred(list_score):
    list_phred = ''
    for score in list_score:
        phred = int(round(score * 93, 0))
        list_phred += chr(phred + 33)
    return list_phred

def run_consensus(path_to_reads):


    if path_to_reads[-3:] == 'sam':
        sequence = read_sam(path_to_reads)
    elif path_to_reads[-5:] == 'fastq':
        sequence = read_fastq(path_to_reads)
    elif path_to_reads[-5:] == 'fasta':
        sequence = read_fasta(path_to_reads)
    write_seq_in_file_with_length('myreads.fasta', sequence, 0, 0)

    # подсчёт средней длины последовательности, чтобы избавиться от слишком коротких и слишком длинных последовательностей
    alligment = read_seq_from_file('myreads.fasta')
    avg_length = sum(len(word) for word in alligment) / len(alligment)
    list_allig = []
    list_need_seq = []
    for j, i in enumerate(alligment):
        if len(i) >= avg_length / 1.6 or len(i) <= avg_length / 0.66:
            list_allig.append(i)
            list_need_seq.append(j)
    write_seq_in_file_with_length('myreads.fasta', list_allig, 0, 0)

    # Стартовые позиции последовательности, у нас всегда 0
    list_startpos = []
    for i in list_allig:
        list_startpos.append(0)

    # Узнавание какая последовательность самая длинная и её индекс
    sequence = read_seq_from_file('myreads.fasta')
    maxi_index = max([len(word) for word in sequence])
    index_max = 0
    for j, i in enumerate(sequence):
        if len(i) == maxi_index:
            index_max = j

    # Основная часть, в которой сначала проверяется подходит ли нам последовательность учитывая стартовую позицию
    # Далее отрезается часть последовательности от indi до indi + 80, где indi это индекс, который увеличивается на половину длины выравнивания 
    # Далее считается количество пропусков и нуклеотидов и заменяется в изначальных последовательностях куски нуклеотидов на выравненные 
    list_already_use = []
    list_lish_gap = []
    indi = 0
    count_nucl = 0
    for start in tqdm(range(0, 20000, 40)):
        list_use_seq = []
        list_use_seq_index = []
        

        for j, seq in enumerate(sequence):
            if count_nucl >= list_startpos[j]:
                if j not in list_already_use:
                    sequence[j] = indi * '-' + seq + (10000 - len(seq) - list_startpos[j]) * '-'
                    list_already_use.append(j)
                gap_flag = False
                for i in seq[indi:indi + 80]:
                    if i != '-':
                        gap_flag = True
                        continue
                if gap_flag:
                    list_use_seq.append(sequence[j])
                    list_use_seq_index.append(j)    

        write_seq_in_file_with_length_and_name(f'allig/reads_{indi}-{indi + 80}.fasta', list_use_seq, list_use_seq_index, indi, indi + 80)
        muscle(f'allig/reads_{indi}-{indi + 80}.fasta', f'allig/allig_{indi}-{indi + 80}.fasta')

        list_seq_with_insert = []
        sorted_seqs = sort_fasta(f'allig/allig_{indi}-{indi + 80}.fasta')
        for rec in sorted_seqs:
            list_seq_with_insert.append(str(rec.seq))
        write_seq_in_file_with_length_and_name(f'allig/allig_{indi}-{indi + 80}.fasta', list_seq_with_insert, list_use_seq_index,  0, 0)

        alligment = read_seq_from_file(f'allig/allig_{indi}-{indi + 80}.fasta')



        tmp = int(0.5 * len(alligment[0]))

        list_allig_index = []
        list_count_gap = []
        for j, seq in enumerate(alligment):
            tmp2 = seq[:tmp].count('-')
            list_count_gap.append(tmp2)
            list_allig_index.append(tmp - tmp2)
        
        for i in range(len(alligment)):
            alligment[i] = alligment[i][:tmp]

        for j, i in enumerate(list_use_seq_index):
            if i == index_max:
                count_nucl += int(list_allig_index[j] * 0.93)
            sequence[i] = sequence[i][:indi] + alligment[j] + sequence[i][list_allig_index[j] + indi:]
        indi += tmp

    # Дополнение строк до одинаковой длины
    list_length = []
    for i in sequence:
        list_length.append(len(i))
    maxi = max(list_length)
    for j, i in enumerate(sequence):
        while len(sequence[j]) < maxi:
            sequence[j] = sequence[j] + '-'

    write_seq_in_file_with_length('allig_seq.fasta', sequence, 0, 0)

    # Постройка консенсуса 
    list_consensus = []
    alligment = read_seq_from_file('allig_seq.fasta')
    count_nucl = per_nucl_in_cal(alligment)
    count_gap = count_gaps (alligment)
    alligment = delete_gap(alligment, count_gap, count_nucl, 50 / 100)
    # Запись множественного выравнивания в файл
    write_seq_in_file_with_length('multiple_alignment.fasta', alligment, 0, 0)
    pi = start_probability (alligment)
    a = matrix_of_MDI(alligment)
    A = probability_of_trans(alligment, a)
    count = count_nucl_in_cal(alligment)
    b = percetn_of_nucl_in_cal(count)
    consensuss = consensus(alligment, b)
    # Запись консенсуса в файл
    write_seq_in_file_with_length('consensus.fasta', [consensuss], 0, 0)
    # Расчёт качества
    consi = read_seq_from_file('consensus.fasta')
    multiple_allig = read_seq_from_file('multiple_alignment.fasta')
    list_seq = consi + multiple_allig
    write_seq_in_file_with_length('consensus_and_multiple.fasta', list_seq, 0, 0)
    list_prob_quality = quality(list_seq)
    phred_string = convert_to_phred(list_prob_quality)
    write_fastq('consensus.fastq', consi[0][:-100], phred_string)
    return