import subprocess
from statistics import median
from Bio import SeqIO
import subprocess


def read_fastq(path_to_reads):
    """
    Читает последовательности из файла FASTQ.
    Аргументы:
        path_to_reads (str): Путь к файлу FASTQ.
    Возвращает:
        list: Список последовательностей.
    """
    sequences = []
    with open(path_to_reads, "r") as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            sequence = lines[i+1].strip()
            sequences.append(sequence)
    return sequences

def write_fastq(file_path, sequence, quality):
    """
    Записывает последовательность и качество в файл FASTQ.
    Аргументы:
        file_path (str): Путь к выходному файлу.
        sequence (str): Последовательность нуклеотидов.
        quality (str): Качество для последовательности.
    """
    with open(file_path, "w") as file:
        file.write("@consensus\n")
        file.write(sequence + "\n")
        file.write("+\n")
        file.write(quality + "\n")

def fastq_to_fasta(path_fastq, path_fasta):
    """
    Преобразует файл FASTQ в формат FASTA.
    Аргументы:
        path_fastq (str): Путь к входному файлу FASTQ.
        path_fasta (str): Путь к выходному файлу FASTA.
    """
    with open(path_fasta, "w") as output:
        sequences = SeqIO.parse(path_fastq, "fastq")
        SeqIO.write(sequences, output, "fasta")

def read_seq_from_file(path):
    """
    Читает последовательности из файла FASTA.
    Аргументы:
        path (str): Путь к файлу FASTA.
    Возвращает:
        list: Список последовательностей.
    """
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
        if seq:
            sequences.append(seq)
    return sequences

def write_seq_in_file_with_length(name_file, list_seq, start, end):
    """
    Записывает срезы последовательностей в файл.
    Аргументы:
        name_file (str): Имя выходного файла.
        list_seq (list): Список последовательностей.
        start (int): Начало среза.
        end (int): Конец среза (не включительно).
    """
    with open(name_file, "w") as f:
        for j, i in enumerate(list_seq):
            header = f">sequence_id{j}\n"
            f.write(header)
            if start == end:
                f.write(i + '\n')
            else:
                f.write(i[start:end] + '\n')

def muscle(path_in, path_out, muscle_bin_full_path):
    """
    Запускает программу MUSCLE с штрафом за пропуск равным 2.
    Аргументы:
        path_in (str): Путь к входному файлу.
        path_out (str): Путь к выходному файлу.
        muscle_bin_full_path (str): Путь к исполняемому файлу MUSCLE.
    """
    command = f"{muscle_bin_full_path} -gapopen -2 -in {path_in} -out {path_out} -quiet"
    subprocess.run(command, shell=True)

def muscle_with_gap1000(path_in, path_out, muscle_bin_full_path):
    """
    Запускает программу MUSCLE с штрафом за пропуск равным 1000.
    Аргументы:
        path_in (str): Путь к входному файлу.
        path_out (str): Путь к выходному файлу.
        muscle_bin_full_path (str): Путь к исполняемому файлу MUSCLE.
    """
    command = f"{muscle_bin_full_path} -gapopen -1000 -in {path_in} -out {path_out} -quiet"
    subprocess.run(command, shell=True)

def read_fasta_start_pos(path_to_reads):
    """
    Читает последовательности из файла FASTA.
    Аргументы:
        path_to_reads (str): Путь к файлу FASTA.
    Возвращает:
        list: Список последовательностей.
    """
    sequences = []
    with open(path_to_reads, "r") as file:
        sequence_lines = []
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_lines:
                    sequences.append(''.join(sequence_lines))
                    sequence_lines = []
            else:
                sequence_lines.append(line)
        if sequence_lines:
            sequences.append(''.join(sequence_lines))
    return sequences

def complement(sequence):
    """
    Возвращает обратный комплемент ДНК-последовательности.
    Аргументы:
        sequence (str): ДНК-последовательность.
    Возвращает:
        str: Обратный комплемент последовательности.
    """
    comp_bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = ''.join(comp_bases[base] for base in sequence)
    return seq[::-1]

def read_sam(path_to_reads):
    """
    Читает последовательности из файла SAM.
    Аргументы:
        path_to_reads (str): Путь к файлу SAM.
    Возвращает:
        list: Список последовательностей.
    """
    sequence = []
    with open(path_to_reads, "r") as sam_file:
        for line in sam_file:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            sequence.append(parts[9])
    return sequence

def min_length_string(strings):
    """
    Возвращает строку минимальной длины из списка строк.
    Аргументы:
        strings (list): Список строк.
    Возвращает:
        str: Строка минимальной длины.
    """
    if not strings:
        return None
    return min(strings, key=len)

def extract_substring(string, left_offset, right_offset, shift=0):
    """
    Возвращает подстроку из строки с учетом отступов от центра и смещения.
    Аргументы:
        string (str): Исходная строка.
        left_offset (int): Отступ влево от центра.
        right_offset (int): Отступ вправо от центра.
        shift (int): Смещение относительно центра (по умолчанию 0).
    Возвращает:
        str: Подстрока.
    """
    length = len(string)
    center = length // 2 + shift
    start = max(0, center - left_offset)
    end = min(length, center + right_offset)
    return string[start:end]



def run_minimap228(long_sequence_file, short_sequence_file, output_file, path_to_minimap2):
    """
    Запускает minimap2 версии 2.28 для выравнивания последовательностей.
    Аргументы:
        long_sequence_file (str): Путь к файлу с длинной последовательностью.
        short_sequence_file (str): Путь к файлу с короткой последовательностью.
        output_file (str): Путь к выходному файлу.
    """
    command = [
        path_to_minimap2,
        "--secondary=no",
        "-t", "20",
        "-k", "7",
        "-w", "1",
        "-F", "5000",
        long_sequence_file,
        short_sequence_file
    ]
    with open(output_file, "w") as outfile:
        subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, text=True)

def parse_paf(file_path):
    """
    Парсит PAF файл и извлекает информацию о выравниваниях.
    Аргументы:
        file_path (str): Путь к PAF файлу.
    Возвращает:
        list: Список словарей с информацией о выравниваниях.
    """
    alignments = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            query_name = fields[0]
            length_1 = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            positive_or_negative = fields[4]
            target_name = fields[5]
            length_2 = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            matching_bases = int(fields[9])
            alignment_length = int(fields[10])
            mapq = int(fields[11])
            alignment_info = {
                'query_name': query_name,
                'length_1': length_1,
                'length_2': length_2,
                'query_start': query_start,
                'query_end': query_end,
                'positive_or_negative': positive_or_negative,
                'target_name': target_name,
                'target_start': target_start,
                'target_end': target_end,
                'matching_bases': matching_bases,
                'alignment_length': alignment_length,
                'mapq': mapq
            }
            alignments.append(alignment_info)
    return alignments

def compute_shifts(alignment_data):
    """
    Вычисляет смещения между длинными и короткими последовательностями.
    Аргументы:
        alignment_data (list): Список словарей с данными выравнивания.
    Возвращает:
        list: Список словарей с номером последовательности и ее смещением.
    """
    shifts = []
    for idx, alignment in enumerate(alignment_data):
        shift = alignment['start_short'] - alignment['start_long']
        shifts.append({
            'sequence_number': idx + 1,
            'shift': shift
        })
    return shifts

def adjust_shifts(shifts):
    """
    Приводит все смещения к положительным значениям.
    Аргументы:
        shifts (list): Список словарей с номерами последовательностей и их смещениями.
    Возвращает:
        list: Список словарей с обновленными смещениями.
    """
    min_shift = min(item['shift'] for item in shifts)
    adjustment = -min_shift
    adjusted_shifts = []
    for item in shifts:
        new_shift = item['shift'] + adjustment
        adjusted_shifts.append({
            'sequence_number': item['sequence_number'],
            'shift': new_shift
        })
    return adjusted_shifts

def compute_median_shifts(data):
    """
    Вычисляет медианные смещения для каждой последовательности.
    Аргументы:
        data (list): Список списков со смещениями.
    Возвращает:
        dict: Словарь с медианными смещениями по номеру последовательности.
    """
    shift_values = {}
    for sequence in data:
        for item in sequence:
            seq_num = item['sequence_number']
            shift = item['shift']
            if seq_num not in shift_values:
                shift_values[seq_num] = []
            shift_values[seq_num].append(shift)
    median_shifts = {}
    for seq_num, shifts in shift_values.items():
        median_shifts[seq_num] = int(median(shifts))
    return median_shifts

def transpose(matrix):
    """
    Транспонирует матрицу.
    Аргументы:
        matrix (list of lists): Исходная матрица.
    Возвращает:
        list of lists: Транспонированная матрица.
    """
    return [list(row) for row in zip(*matrix)]

def parse_startpos(qname):
    """
    Извлекает значение startpos из поля QNAME в SAM файле.
    Аргументы:
        qname (str): Значение поля QNAME.
    Возвращает:
        str: Значение startpos или None, если не найдено.
    """
    qname_fields = qname.split(';')
    for field in qname_fields:
        field = field.strip()
        if field.startswith('startpos='):
            startpos_value = field[len('startpos='):]
            return startpos_value
    return None

def extract_startpos_from_sam(sam_file_path):
    """
    Извлекает значения startpos из SAM файла.
    Аргументы:
        sam_file_path (str): Путь к SAM файлу.
    Возвращает:
        list: Список значений startpos.
    """
    list_startpos = []
    with open(sam_file_path, 'r') as sam_file:
        for line in sam_file:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            qname = fields[0]
            startpos = parse_startpos(qname)
            if startpos is not None:
                list_startpos.append(startpos)
    return list_startpos

def determine_orientations(signs):
    """
    Определяет ориентации последовательностей на основе знаков выравниваний.
    Аргументы:
        signs (list of lists): Матрица знаков выравниваний.
    Возвращает:
        tuple: Списки номеров прямых и обратных последовательностей.
    """
    sign_map = {'+': 1, '-': -1}
    N = len(signs)
    signs_int = []
    for row in signs:
        signs_int.append([sign_map[c] for c in row])

    max_agreements = -1
    best_labels = None

    for assignment in range(1 << N):
        labels = [1 if (assignment >> i) & 1 else -1 for i in range(N)]
        agreements = 0
        for i in range(N):
            for j in range(N):
                if signs_int[i][j] == labels[i] * labels[j]:
                    agreements += 1
        if agreements > max_agreements:
            max_agreements = agreements
            best_labels = labels

    direct_indices = [idx + 1 for idx, label in enumerate(best_labels) if label == 1]
    reverse_indices = [idx + 1 for idx, label in enumerate(best_labels) if label == -1]

    return direct_indices, reverse_indices

def start_pos(file_reads, file_out, path_to_minimap2, path_to_outdir):
    """
    Выравнивает все последовательности из файла file_reads и сохраняет результат в file_out.
    Аргументы:
        file_reads (str): Путь к файлу с исходными последовательностями в формате FASTA.
        file_out (str): Путь к файлу для сохранения выровненных последовательностей.
    Возвращает:
        tuple: (list_med_real, list_all_positive_or_negative)
    """
    # Выравнивание всей длины наименьшего на все остальные
    reads = read_fasta_start_pos(file_reads)

    list_all_border = []
    list_all_positive_or_negative = []

    for index1, read_align in enumerate(reads):
        list_border = []
        list_positive_or_negative = []
        write_seq_in_file_with_length(f'{path_to_outdir}short_sequence.fasta', [read_align], 0, 1111111)
        for index, read in enumerate(reads):
            write_seq_in_file_with_length(f'{path_to_outdir}long_sequence.fasta', [read], 0, 1111111)
            run_minimap228(long_sequence_file=f"{path_to_outdir}long_sequence.fasta", 
                           short_sequence_file=f"{path_to_outdir}short_sequence.fasta", 
                           output_file=f"{path_to_outdir}alignments_{index1}_{index}.paf", 
                           path_to_minimap2=path_to_minimap2)
            paf_file_path = f'{path_to_outdir}alignments_{index1}_{index}.paf'  # Укажите путь к вашему PAF файлу
            alignments = parse_paf(paf_file_path)
            #print(alignments)
            if alignments:
                for alignment in alignments:
                    list_border.append({
                        "start_long": alignment['target_start'], 
                        "end_long": alignment['target_end'],
                        "length_long": alignment['length_2'],
                        "start_short": alignment['query_start'], 
                        "end_short": alignment['query_end'],
                        "length_short": alignment['length_1'],
                    })
                    list_positive_or_negative.append(alignment['positive_or_negative'])
            else:
                continue
        list_all_border.append(list_border)
        list_all_positive_or_negative.append(list_positive_or_negative)

    #for i in list_all_border:
        #print(i)

    #for i in list_all_positive_or_negative:
        #print(i)

    list_shifts = []
    # Вычисляем смещения
    for i in list_all_border:
        list_shifts.append(compute_shifts(i))

    #print(list_shifts)

    # Применяем функцию
    list_adjust_shifts = []
    for i in list_shifts:
        list_adjust_shifts.append(adjust_shifts(i))

    #print(list_adjust_shifts)

    #for i in list_adjust_shifts:
    #   for j in i:
            #print(j['shift'], sep=' ', end=' ')
        #print()

    # Вызов функции и вывод результатов
    median_shifts = compute_median_shifts(list_adjust_shifts)
    list_med_real = []
    for index, seq_num in enumerate(sorted(median_shifts.keys())):
        #print(f"Sequence Number {seq_num}: Median Shift = {median_shifts[seq_num]}")
        list_med_real.append([int(median_shifts[seq_num])])

    sequence = read_fasta_start_pos(file_reads)

    for i in range(len(sequence)):
        sequence[i] = '-' * median_shifts[i + 1] + sequence[i]

    write_seq_in_file_with_length(file_out, sequence, 0, 11111111)

    #for i in sequence:
        #print(i)

    return list_med_real, list_all_positive_or_negative


def filter_fasta_by_length(input_file, min_length=50):
    """
    Фильтрует последовательности в FASTA-файле по минимальной длине.

    Аргументы:
        input_file (str): Путь к входному FASTA-файлу.
        min_length (int): Минимальная допустимая длина последовательности.
    """
    sequences = []
    with open(input_file, 'r') as infile:
        header = ''
        seq = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header and len(seq) >= min_length:
                    sequences.append((header, seq))
                header = line
                seq = ''
            else:
                seq += line
        if header and len(seq) >= min_length:
            sequences.append((header, seq))
    with open(input_file, 'w') as outfile:
        for header, seq in sequences:
            outfile.write(f'{header}\n')
            for i in range(0, len(seq), 60):
                outfile.write(f'{seq[i:i+60]}\n')
