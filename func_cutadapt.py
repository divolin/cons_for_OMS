import subprocess
import os

def read_fasta(path_to_reads):
    """
    Читает последовательности из файла FASTA.
    Аргументы:
        path_to_reads (str): Путь к файлу FASTA.
    Возвращает:
        list: Список последовательностей.
    """
    sequence = []
    with open(path_to_reads, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                sequence.append(line)
    return sequence

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


def complement(sequence):
    """
    Возвращает обратную ДНК-последовательность.
    Аргументы:
        sequence (str): ДНК-последовательность.
    Возвращает:
        str: Обратный комплемент последовательности.
    """
    comp_bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = ''.join(comp_bases[base] for base in sequence)
    return seq[::-1]

def write_seq_in_file_with_length(name_file, list_seq, start, end):
    """
    Записывает срезы последовательностей в файл.
    Аргументы:
        name_file (str): Имя файла для записи.
        list_seq (list): Список последовательностей.
        start (int): Начало среза.
        end (int): Конец среза (не включительно).
    """
    with open(name_file, "w") as f:
        for j, i in enumerate(list_seq):
            header = f">sequence_id{j}\n"
            f.write(header)
            if start == end:
                f.write(i +'\n')
            else:
                f.write(i[start:end] +'\n')

def parse_cutadapt_summary(file_path):
    """
    Парсит отчет cutadapt для извлечения обрезанных длин слева и справа.
    Аргументы:
        file_path (str): Путь к файлу отчета cutadapt.
    Возвращает:
        tuple: (left_trimmed, right_trimmed)
    """
    left_trimmed = None
    right_trimmed = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if 'Overview of removed sequences at 5\' end' in line:
                for j in range(i + 2, len(lines)):
                    if lines[j].strip() == '':
                        break
                    length, count, _, max_err, _ = lines[j].split('\t')
                    left_trimmed = int(length)
            if 'Overview of removed sequences at 3\' end' in line:
                for j in range(i + 2, len(lines)):
                    if lines[j].strip() == '':
                        break
                    length, count, _, max_err, _ = lines[j].split('\t')
                    right_trimmed = int(length)
    return left_trimmed, right_trimmed

def parse_cutadapt_summary_one(file_path):
    """
    Парсит отчет cutadapt для извлечения обрезанной длины слева.
    Аргументы:
        file_path (str): Путь к файлу отчета cutadapt.
    Возвращает:
        int: Обрезанная длина слева.
    """
    left_trimmed = None
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if 'Overview of removed sequences' in line:
                for j in range(i + 2, len(lines)):
                    if lines[j].strip() == '':
                        break
                    length, count, _, max_err, _ = lines[j].split('\t')
                    left_trimmed = int(length)
    return left_trimmed

def read_sam(path_to_reads):
    """
    Читает последовательности из SAM-файла.
    Аргументы:
        path_to_reads (str): Путь к SAM-файлу.
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








def process_reads_with_cutadapt(long_read_file, adapter_start, adapter_end, path_to_outdir, g_a, left_right):
    """
    Обрабатывает файл длинных прочтений с помощью инструмента cutadapt, выполняя итеративное обрезание адаптеров
    и обновление последовательностей.

    Аргументы:
        long_read_file (str): Путь к входному файлу с длинными прочтениями в формате FASTA.
        adapter_start (str): Последовательность адаптера в начале прочтения.
        adapter_end (str): Последовательность адаптера в конце прочтения.
        path_to_outdir (str): Путь к выходной директории для сохранения промежуточных файлов.

    Функция использует внешние инструменты и функции:
        - 'cutadapt': инструмент для обрезания адаптеров.
        - 'parse_cutadapt_summary': функция для парсинга вывода cutadapt (должна быть реализована).
        - 'read_fasta': функция для чтения FASTA файлов (должна быть реализована).
        - 'write_seq_in_file_with_length': функция для записи последовательностей в файл (должна быть реализована).

    Возвращает:
        None
    """
    
    str_cutadapt = f'{adapter_start}...{adapter_end}'

    count = 0  
    flag = True  

    while flag:

        cmd = [
            'cutadapt',
            g_a, str_cutadapt,  
            '-e', '0.3',  
            '-o', os.path.join(path_to_outdir, 'output.fasta'),  
            long_read_file  
        ]

        with open(os.path.join(path_to_outdir, 'output.txt'), 'w') as outfile:
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.DEVNULL)

        left_trimmed, right_trimmed = parse_cutadapt_summary(os.path.join(path_to_outdir, 'output.txt'))

        if left_trimmed:
            if right_trimmed:
                sequence_after_cutadapt = read_fasta(os.path.join(path_to_outdir, 'output.fasta'))
                if sequence_after_cutadapt:
                    sequence_after_cutadapt = sequence_after_cutadapt[0]
                    if len(sequence_after_cutadapt) < 10:
                        count += 1
                        if count == 4:
                            flag = False
                            break
                    else:
                        if left_right == 'yes':
                            left_trimmed -= len(adapter_start)
                            if left_trimmed < 0:
                                left_trimmed = 0
                            right_trimmed += len(adapter_start)
                        else:
                            right_trimmed += len(adapter_start)
                else:
                    count += 1
                    if count == 4:
                        flag = False
                        break
            else:
                flag = False
                break
        else:
            flag = False
            break

        sequence_after_cutadapt = read_fasta(os.path.join(path_to_outdir, 'output.fasta'))[0]

        all_seq_after_cutadapt = read_fasta(os.path.join(path_to_outdir, 'reads.fasta'))
        all_seq_after_cutadapt.append(sequence_after_cutadapt)

        write_seq_in_file_with_length(os.path.join(path_to_outdir, 'reads.fasta'), all_seq_after_cutadapt, 0, 0)

        long_sequence = read_fasta(long_read_file)[0]

        long_sequence = (
            long_sequence[0:left_trimmed] +
            long_sequence[len(sequence_after_cutadapt) + left_trimmed + len(adapter_start):]
        )

        write_seq_in_file_with_length(long_read_file, [long_sequence], 0, 0)




def process_reads_with_cutadapt_one(long_read_file, adapter_start, adapter_end, path_to_outdir, g_a, shift):
    """
    Обрабатывает файл длинных прочтений с помощью инструмента cutadapt, выполняя итеративное обрезание адаптеров
    и обновление последовательностей.

    Аргументы:
        long_read_file (str): Путь к входному файлу с длинными прочтениями в формате FASTA.
        adapter_start (str): Последовательность адаптера в начале прочтения.
        adapter_end (str): Последовательность адаптера в конце прочтения.
        path_to_outdir (str): Путь к выходной директории для сохранения промежуточных файлов.

    Функция использует внешние инструменты и функции:
        - 'cutadapt': инструмент для обрезания адаптеров.
        - 'parse_cutadapt_summary': функция для парсинга вывода cutadapt (должна быть реализована).
        - 'read_fasta': функция для чтения FASTA файлов (должна быть реализована).
        - 'write_seq_in_file_with_length': функция для записи последовательностей в файл (должна быть реализована).

    Возвращает:
        None
    """

    if g_a == '-g':
        str_cutadapt = f'{adapter_end};rightmost'
    else:
        str_cutadapt = adapter_end

    count = 0  
    flag = True  

    while flag:
        cmd = [
            'cutadapt',
            g_a, str_cutadapt,  
            '-e', '0.3', 
            '-o', os.path.join(path_to_outdir, 'output.fasta'),  
            long_read_file  
        ]

        with open(os.path.join(path_to_outdir, 'output.txt'), 'w') as outfile:
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.DEVNULL)

        left_trimmed = parse_cutadapt_summary_one(os.path.join(path_to_outdir, 'output.txt'))

        if left_trimmed:
            sequence_after_cutadapt = read_fasta(os.path.join(path_to_outdir, 'output.fasta'))
            if sequence_after_cutadapt:
                sequence_after_cutadapt = sequence_after_cutadapt[0]
                if len(sequence_after_cutadapt) < 10:
                    count += 1
                    if count == 4:
                        flag = False
                        break
                all_seq_after_cutadapt = read_fasta(os.path.join(path_to_outdir, 'reads.fasta'))
                all_seq_after_cutadapt.append(sequence_after_cutadapt)

                write_seq_in_file_with_length(os.path.join(path_to_outdir, 'reads.fasta'), all_seq_after_cutadapt, 0, 0)

                long_sequence = read_fasta(long_read_file)[0]

                if shift == 0:
                    long_sequence = (
                        long_sequence[0:left_trimmed] +
                        long_sequence[len(sequence_after_cutadapt) + left_trimmed:]
                    )
                else:
                    left_board = left_trimmed - count * int(len(adapter_end) * 0.9)
                    if left_board < 0:
                        left_board = 0
                    right_board = len(sequence_after_cutadapt) + left_trimmed
                    long_sequence = (
                        long_sequence[0:left_board] +
                        long_sequence[right_board:]
                    )

                write_seq_in_file_with_length(long_read_file, [long_sequence], 0, 0)
            else:
                count += 1
                if count == 4:
                    flag = False
                    break
        else:
            flag = False
            break
        






def cutadapt_func(adapter_start, adapter_end, long_read_file, path_to_outdir, output_folder):
    """
    Выполняет обрезку адаптеров с помощью cutadapt на файле длинных чтений.
    Аргументы:
        adapter_start (str): Начальный адаптер.
        adapter_end (str): Конечный адаптер.
        long_read_file (str): Путь к файлу длинных чтений в формате FASTA.
    """
    with open(os.path.join(output_folder, 'reads.fasta'), 'w') as file:
        pass
    
    process_reads_with_cutadapt(
        long_read_file=long_read_file,
        adapter_start=adapter_start,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-a',
        left_right='yes')

    process_reads_with_cutadapt(
        long_read_file=long_read_file,
        adapter_start=adapter_end,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-a',
        left_right='no')
    
    process_reads_with_cutadapt_one(
        long_read_file=long_read_file,
        adapter_start=adapter_end,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-g',
        shift=0)
        
    process_reads_with_cutadapt(
        long_read_file=long_read_file,
        adapter_start=adapter_start,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-a',
        left_right='no')
    
    process_reads_with_cutadapt(
        long_read_file=long_read_file,
        adapter_start=adapter_start,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-a',
        left_right='no')
                
    process_reads_with_cutadapt_one(
        long_read_file=long_read_file,
        adapter_start=adapter_end,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-a',
        shift=0)
    
    process_reads_with_cutadapt_one(
        long_read_file=long_read_file,
        adapter_start=adapter_end,
        adapter_end=adapter_end,
        path_to_outdir=output_folder,
        g_a='-g',
        shift=10)

    long_sequence = read_fasta(f'{long_read_file}')[0]
    all_seq_after_cutadapt = read_fasta(os.path.join(output_folder, 'reads.fasta'))
    all_seq_after_cutadapt.append(long_sequence)
    write_seq_in_file_with_length(os.path.join(output_folder, 'reads.fasta'), all_seq_after_cutadapt, 0, 111111111)
