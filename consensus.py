from func_cutadapt import *
from func_start_pos import *
from func_for_BioUML import *
import traceback
import shutil
import os
import sys
import argparse
from pathlib import Path
import time

def delete_file(file):
    try:
        os.remove(file)
    except:
        pass

def remove_initial_gaps(file_path):
    """
    Удаляет ведущие гэпы '-' во всех последовательностях, 
    если во всех строках на данной позиции стоит '-'.

    :param sequences: список строк (ДНК-последовательностей)
    :return: список строк без начальных общих гэпов
    """
    sequences = read_fasta_start_pos(file_path)

    if not sequences or all(len(seq) == 0 for seq in sequences):
        return sequences

    min_len = min(len(seq) for seq in sequences)

    col_index = 0

    while col_index < min_len:
        if all(seq[col_index] == '-' for seq in sequences):
            col_index += 1
        else:
            break
    trimmed_sequences = [seq[col_index:] for seq in sequences]
    write_seq_in_file_with_length(file_path, trimmed_sequences, 0, 0)
    return



def build_consensus(adapter_fasta, 
                    path_to_reads, 
                    path_to_outdir, 
                    muscle_bin_full_path, 
                    name_consensus,
                    path_to_minimap2):
    start_time = time.time()
    output_folder = Path(path_to_outdir)
    output_folder.mkdir(parents=True, exist_ok=True)

    output_folder = Path(os.path.join(path_to_outdir, name_consensus.removesuffix('.fastq')))
    output_folder.mkdir(parents=True, exist_ok=True)

    print(f'Обработка файла: {path_to_reads}')
    with open(os.path.join(output_folder, "log.txt"), "a") as f:
        print('Идёт обрезка адаптеров', file=f)
    try:
        adapter_start = read_fasta(adapter_fasta)[0]
        adapter_end = read_fasta(adapter_fasta)[0]

        if path_to_reads.endswith('.sam'):
            sequence = read_sam(path_to_reads)
        elif path_to_reads.endswith('.fastq'):
            sequence = read_fastq(path_to_reads)
        elif path_to_reads.endswith('.fasta'):
            sequence = read_fasta_start_pos(path_to_reads)
        else:
            with open(os.path.join(output_folder, "log.txt"), "a") as f:
                print(f'Пропуск файла {path_to_reads}: неподдерживаемое расширение.', file=f)
            return  

        read_length_poly = len(sequence[0])

        long_read_file = os.path.join(output_folder, 'long_read.fasta')
        write_seq_in_file_with_length(long_read_file, sequence, 0, 11111111111)

        

        cutadapt_func(adapter_start=adapter_start,
                    adapter_end=adapter_end,
                    long_read_file=long_read_file,
                    path_to_outdir=path_to_outdir,
                    output_folder=output_folder)
        delete_file(os.path.join(output_folder, 'output.fasta'))
        delete_file(os.path.join(output_folder, 'output.txt'))

        input_fasta = os.path.join(output_folder, 'reads.fasta')  
        filter_fasta_by_length(input_fasta)

        reads_after_cutadapt = read_fasta_start_pos(os.path.join(output_folder, 'reads.fasta')  )
        count_reads = len(reads_after_cutadapt)
        
        average_length = round(sum(len(s) for s in reads_after_cutadapt) / count_reads)
        obr_time = time.time()
        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            print(f'Обрезка адаптеров окончена, получилось: {count_reads} прочтений со средней длиной: {average_length}', file=f)
            
            print(f"Время обрезки адаптеров {obr_time - start_time}", file=f)
    except:
        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            traceback.print_exc(file=f)
        
            print('Ошибка в обрезке адаптеров', file=f)
            print('--------------------------------------------------', file=f)
        return 0

    with open(os.path.join(output_folder, "log.txt"), "a") as f:
        print('Идёт поиск стартовых позиций', file=f)
    try:
        alligment = read_seq_from_file(os.path.join(output_folder, 'reads.fasta'))
        avg_length = sum(len(word) for word in alligment) / len(alligment)
        lower_bound = avg_length * 0.7
        upper_bound = avg_length * 1.3
        list_allig = [seq for seq in alligment if lower_bound <= len(seq) <= upper_bound]
        write_seq_in_file_with_length(os.path.join(output_folder, 'reads.fasta'), list_allig, 0, 0)

        
        list_all_positive_or_negative = start_pos_positive_or_negative(file_reads=os.path.join(output_folder, 'reads.fasta'), 
                                                                 path_to_minimap2=path_to_minimap2, 
                                                                 path_to_outdir=output_folder)
        list_napravlenie_reads = determine_orientations(list_all_positive_or_negative)
        list_all_reads = read_fasta_start_pos(os.path.join(output_folder, 'reads.fasta')  )
        list_reads = [list_all_reads[i - 1] for i in list_napravlenie_reads[0]]
        list_compl_reads = [list_all_reads[i - 1] for i in list_napravlenie_reads[1]]
        write_seq_in_file_with_length(os.path.join(output_folder, 'real_reads.fasta'), list_reads, 0, 0)
        write_seq_in_file_with_length(os.path.join(output_folder, 'real_compl_reads.fasta'), list_compl_reads, 0, 0)


        start_pos(file_reads=os.path.join(output_folder, 'real_reads.fasta'), 
                file_out=os.path.join(output_folder, 'out.fasta'), 
                path_to_minimap2=path_to_minimap2,
                path_to_outdir=output_folder)
        start_pos(file_reads=os.path.join(output_folder, 'real_compl_reads.fasta'), 
                file_out=os.path.join(output_folder, 'out_compl.fasta'), 
                path_to_minimap2=path_to_minimap2,
                path_to_outdir=output_folder)
        
        remove_initial_gaps(os.path.join(output_folder, 'out.fasta'))
        remove_initial_gaps(os.path.join(output_folder, 'out_compl.fasta'))

        for i in range(100):
            for j in range(100):
                delete_file(os.path.join(output_folder, f'alignments_{i}_{j}.paf'))
        st_time = time.time()
        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            print('Поиск стартовых позиций окончен', file=f)
            print(f"Время поиска стартовых позиций {st_time - obr_time}", file=f)
    except Exception as e:
        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            traceback.print_exc(file=f)
            print('Ошибка в поиске стартовых позиций', file=f)
            print('--------------------------------------------------', file=f)
        return 0
    
    with open(os.path.join(output_folder, "log.txt"), "a") as f:
        print('Идёт составление консенсуса', file=f)
    try:
        # Запуск для прямых прочтений
        run_consensus(path_to_reads=os.path.join(output_folder, 'out.fasta'), 
                      path_to_outdir=output_folder,
                      muscle_bin_full_path=muscle_bin_full_path)

        # Запуск для обратных прочтений
        run_consensus_compl(path_to_reads=os.path.join(output_folder, 'out_compl.fasta'), 
                            path_to_outdir=output_folder, 
                            muscle_bin_full_path=muscle_bin_full_path)
    # Получение одного консенсуса из прямого и обратного консенсуса
        consensus_final(
            muscle_bin_full_path=muscle_bin_full_path,
            path_out=output_folder,
            name_consensus=name_consensus
        )
    


        temporary_files = [
            'long_read.fasta', 'long_sequence.fasta', 'short_sequence.fasta', 'predict_1.fasta',
            'tmp.fasta', 'tmp_allig.fasta', 'consensus_compl.fasta'
        ]

        for folder_name in ['allig', 'allig_compl']:
            folder_path = os.path.join(output_folder, folder_name)
            try:
                shutil.rmtree(folder_path)
            except FileNotFoundError:
                pass
            except Exception as e:
                with open(os.path.join(output_folder, "log.txt"), "a") as f:
                    print(f"Произошла ошибка при удалении папки '{folder_path}': {e}", file=f)

        for temp_file in temporary_files:
            delete_file(os.path.join(output_folder, temp_file))

        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            print(f'Составление консенсуса окончено, итоговый консенсус записан по пути {output_folder} в файл {name_consensus}', file=f)

        sequence_path = os.path.join(output_folder, name_consensus)
        sequence = read_fastq(sequence_path)[0]
        end_time = time.time()
        con_time = time.time()
        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            print(f"Изначальная длина полимеразного прочтения {read_length_poly}", file=f)
            print(f"Адаптер начала {adapter_start}", file=f)
            print(f"Адаптер конца {adapter_end}", file=f)
            print(f"После вырезки адаптеров получилось {count_reads} прочтений", file=f)
            print(f"Средняя длина субпрочтений {average_length} нуклеотидов", file=f)
            print(f"Длина консенсуса {len(sequence)} нуклеотидов", file=f)
            print(f"Время построения консенсуса {con_time - st_time}", file=f)
            print(f'Полное время выполнения: {end_time - start_time}', file=f)
            
            print('--------------------------------------------------', file=f)
        
        
        shutil.move(sequence_path, os.path.join(path_to_outdir, name_consensus))
        return 1
    except:
        with open(os.path.join(output_folder, "log.txt"), "a") as f:
            traceback.print_exc(file=f)
            print('Ошибка в составлении консенсуса', file=f)
            print('--------------------------------------------------', file=f)
        return 0








import os
import subprocess

def run_blast(query_fasta, ref_fasta, db_prefix):
    """
    Запускает makeblastdb и blastn для пары (рид, референс) и возвращает
    наибольшее pident (процент совпадений) среди найденных выравниваний.
    Если выравниваний нет, возвращает 0.
    """
    cmd_make_db = f"makeblastdb -in {ref_fasta} -dbtype nucl -out {db_prefix}"
    subprocess.run(cmd_make_db, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    out_file = f"{db_prefix}_blast.out"
    cmd_blastn = (
        f"blastn -query {query_fasta} -db {db_prefix} "
        f"-out {out_file} -outfmt 6 "
    )
    subprocess.run(cmd_blastn, shell=True, check=True)

    best_pident = 0.0
    with open(out_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            pident = float(cols[2])
            bitscore = float(cols[11])
            if pident > 0 and bitscore > 0:
                if bitscore > best_pident:  
                    best_pident = bitscore
                    best_pident_value = pident

    if best_pident > 0:
        return best_pident_value
    else:
        return 0.0

def make_accuracy(reference, reads):
    """
    Вычисляет точность выравнивания последовательности `reads` к референсу `reference`
    с помощью blastn. Возвращает процент совпадений (pident) в лучшем выравнивании.
    Учитываются прямой и обратнокомплементарный варианты референса.
    """

    if reference.split('.')[-1] == 'fasta':
        refe = read_fasta_start_pos(reference)[0]
    elif reference.split('.')[-1] == 'fastq':
        refe = read_fastq(reference)[0]
    else:
        raise ValueError("Формат референса не распознан (нужно .fasta или .fastq)")

    refe_compl = reverse_complement(refe)

    if reads.split('.')[-1] == 'fasta':
        seq = read_fasta_start_pos(reads)[0]
    elif reads.split('.')[-1] == 'fastq':
        seq = read_fastq(reads)[0]
    else:
        raise ValueError("Формат рида не распознан (нужно .fasta или .fastq)")

    write_seq_in_file_with_length('ref.fasta', [refe], 0, 0)
    write_seq_in_file_with_length('ref_compl.fasta', [refe_compl], 0, 0)
    write_seq_in_file_with_length('query.fasta', [seq], 0, 0)

    best_pident_direct = run_blast(
        query_fasta='query.fasta',
        ref_fasta='ref.fasta',
        db_prefix='ref_db'
    )

    best_pident_compl = run_blast(
        query_fasta='query.fasta',
        ref_fasta='ref_compl.fasta',
        db_prefix='ref_compl_db'
    )

    best_pident = max(best_pident_direct, best_pident_compl)

    for ftmp in [
        'ref.fasta', 'ref.nhr', 'ref.nin', 'ref.nsq',
        'ref_db_blast.out', 'ref_db.nhr', 'ref_db.nin', 'ref_db.nsq',
        'ref_compl.fasta', 'ref_compl_db_blast.out',
        'ref_compl_db.nhr', 'ref_compl_db.nin', 'ref_compl_db.nsq',
        'query.fasta'
    ]:
        if os.path.exists(ftmp):
            os.remove(ftmp)

    return round(best_pident, 2)




















if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Обработка последовательностей ДНК.\n* - обязательные параметры')
    parser.add_argument('-a', '--adapter_fasta',
                        required=True,
                        help='*(Обязательный параметр) Путь к файлу адаптера (fasta)')
    parser.add_argument('-r', '--path_to_reads',
                        required=True,
                        help='*(Обязательный параметр) Путь к файлу с прочтениями или к папке с файлами') 
    parser.add_argument('-o', '--path_to_outdir', default=os.path.join(os.getcwd(), 'output'),
                        help='Путь к выходной директории')
    parser.add_argument('-m', '--muscle_bin_full_path', default=os.path.join(os.getcwd(), 'muscle3.8.31_i86linux64'),
                        help='Полный путь к исполняемому файлу MUSCLE')
    parser.add_argument('-s', '--path_to_minimap2', default=os.path.join(os.getcwd(), "minimap2-2.28_x64-linux/minimap2"),
                        help='Полный путь к исполняемому файлу minimap2')

    args = parser.parse_args()

    if os.path.isfile(args.path_to_reads):
        base_name = os.path.splitext(os.path.basename(args.path_to_reads))[0]
        name_consensus = f'{base_name}_consensus.fastq'

        build_consensus(
            adapter_fasta=args.adapter_fasta, 
            path_to_reads=args.path_to_reads, 
            path_to_outdir=args.path_to_outdir, 
            muscle_bin_full_path=args.muscle_bin_full_path, 
            name_consensus=name_consensus,
            path_to_minimap2=args.path_to_minimap2
        )
    elif os.path.isdir(args.path_to_reads):
        supported_extensions = ('.fasta', '.fastq', '.sam')
        files = [os.path.join(args.path_to_reads, f) for f in os.listdir(args.path_to_reads) if f.endswith(supported_extensions)]

        if not files:
            print('В указанной папке нет файлов с расширениями .fasta, .fastq или .sam.')
            sys.exit(1)

        for file_path in files:
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            name_consensus = f'{base_name}_consensus.fastq'

            build_consensus(
                adapter_fasta=args.adapter_fasta, 
                path_to_reads=file_path, 
                path_to_outdir=args.path_to_outdir, 
                muscle_bin_full_path=args.muscle_bin_full_path, 
                name_consensus=name_consensus,
                path_to_minimap2=args.path_to_minimap2
            )
    else:
        print('Параметр -r должен быть файлом или папкой.')
        sys.exit(1)
