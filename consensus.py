from func_cutadapt import *
from func_start_pos import *
from func_for_BioUML import *

import shutil
import os
import sys
import argparse
from pathlib import Path

def delete_file(file):
    try:
        os.remove(file)
    except:
        pass

def build_consensus(adapter_fasta, 
                    path_to_reads, 
                    path_to_outdir, 
                    muscle_bin_full_path, 
                    name_consensus,
                    path_to_minimap2):
    
    # Создаем выходную папку, если она не существует
    output_folder = Path(path_to_outdir)
    output_folder.mkdir(parents=True, exist_ok=True)

    print(f'Обработка файла: {path_to_reads}')
    print('Идёт обрезка адаптеров')
    adapter_start = read_fasta(adapter_fasta)[0]
    adapter_end = read_fasta(adapter_fasta)[0]

    # Определяем тип входного файла и читаем последовательность
    if path_to_reads.endswith('.sam'):
        sequence = read_sam(path_to_reads)
    elif path_to_reads.endswith('.fastq'):
        sequence = read_fastq(path_to_reads)
    elif path_to_reads.endswith('.fasta'):
        sequence = read_fasta(path_to_reads)
    else:
        print(f'Пропуск файла {path_to_reads}: неподдерживаемое расширение.')
        return  # Пропускаем файл, если расширение не поддерживается

    read_length_poly = len(sequence)

    # Создаем временную папку для текущего файла
    temp_dir = os.path.join(path_to_outdir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)


    write_seq_in_file_with_length('long_read.fasta', sequence, 0, 11111111111)

    long_read_file = 'long_read.fasta'

    cutadapt_func(adapter_start=adapter_start,
                  adapter_end=adapter_end,
                  long_read_file=long_read_file,
                  path_to_outdir=path_to_outdir)
    delete_file(f'{path_to_outdir}output.fasta')
    delete_file(f'{path_to_outdir}output.txt')

    input_fasta = f'{path_to_outdir}reads.fasta'  
    filter_fasta_by_length(input_fasta)

    reads_after_cutadapt = read_fasta_start_pos(f'{path_to_outdir}reads.fasta')
    count_reads = len(reads_after_cutadapt)
    average_length = round(sum(len(s) for s in reads_after_cutadapt) / count_reads)

    print(f'Обрезка адаптеров окончена, получилось: {count_reads} прочтений со средней длиной: {average_length}')
    print('Идёт поиск стартовых позиций')

    list_med_real, list_all_positive_or_negative = start_pos(file_reads=f'{path_to_outdir}reads.fasta', 
                                                             file_out=f'{path_to_outdir}predict_1.fasta', 
                                                             path_to_minimap2=path_to_minimap2, 
                                                             path_to_outdir=path_to_outdir)

    list_napravlenie_reads = determine_orientations(list_all_positive_or_negative)
    list_all_reads = read_fasta_start_pos(f'{path_to_outdir}reads.fasta')

    list_reads = [list_all_reads[i - 1] for i in list_napravlenie_reads[0]]
    list_compl_reads = [list_all_reads[i - 1] for i in list_napravlenie_reads[1]]

    write_seq_in_file_with_length(f'{path_to_outdir}real_reads.fasta', list_reads, 0, 1111111)
    write_seq_in_file_with_length(f'{path_to_outdir}real_compl_reads.fasta', list_compl_reads, 0, 1111111)

    start_pos(file_reads=f'{path_to_outdir}real_reads.fasta', 
              file_out=f'{path_to_outdir}out.fasta', 
              path_to_minimap2=path_to_minimap2,
              path_to_outdir=path_to_outdir)
    start_pos(file_reads=f'{path_to_outdir}real_compl_reads.fasta', 
              file_out=f'{path_to_outdir}out_compl.fasta', 
              path_to_minimap2=path_to_minimap2,
              path_to_outdir=path_to_outdir)

    for i in range(50):
        for j in range(50):
            delete_file(f'{path_to_outdir}alignments_{i}_{j}.paf')
    print('Поиск стартовых позиций окончен')

    print('Идёт составление консенсуса')
    # Запуск для прямых прочтений
    run_consensus(path_to_reads=f'{path_to_outdir}out.fasta', 
                  path_to_outdir=path_to_outdir,
                  muscle_bin_full_path=muscle_bin_full_path)

    # Запуск для обратных прочтений
    run_consensus_compl(path_to_reads=f'{path_to_outdir}out_compl.fasta', 
                        path_to_outdir=path_to_outdir, 
                        muscle_bin_full_path=muscle_bin_full_path)
# Получение одного консенсуса из прямого и обратного консенсуса
    consensus_final(
        muscle_bin_full_path=muscle_bin_full_path,
        path_out=path_to_outdir,
        name_consensus=name_consensus
    )

    # Очистка временных файлов
    temporary_files = [
        'long_read.fasta', 'long_sequence.fasta', 'short_sequence.fasta', 'predict_1.fasta',
        'out.fasta', 'out_compl.fasta', 'allig_seq.fasta', 'real_reads.fasta',
        'allig_seq_compl.fasta', 'myreads.fasta', 'tmp.fasta', 'tmp_allig.fasta',
        'real_compl_reads.fasta', 'reads.fasta', 'multiple_alignment_compl.fasta', 'multiple_alignment.fasta',
        'consensus_compl.fastq', 'consensus.fastq'
    ]

    for folder_name in ['allig', 'allig_compl']:
        folder_path = os.path.join(temp_dir, folder_name)
        try:
            shutil.rmtree(folder_path)
        except FileNotFoundError:
            pass
        except Exception as e:
            print(f"Произошла ошибка при удалении папки '{folder_path}': {e}")

    for temp_file in temporary_files:
        delete_file(os.path.join(temp_dir, temp_file))

    # Удаляем временную папку
    shutil.rmtree(temp_dir)

    print(f'Составление консенсуса окончено, итоговый консенсус записан по пути {output_folder} в файл {name_consensus}')

    sequence = read_fastq('/data/output/consensus_final.fastq')[0]

    print(f"Изначальная длина полимеразного прочтения {read_length_poly}")
    print(f"Адаптер начала {adapter_start}")
    print(f"Адаптер конца {adapter_end}")
    print(f"После вырезки адаптеров получилось {count_reads} прочтений")
    print(f"Средняя длина субпрочтений {average_length} нуклеотидов")
    print(f"Длина консенсуса {len(sequence)} нуклеотидов")

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

    # Проверяем, является ли path_to_reads файлом или папкой
    if os.path.isfile(args.path_to_reads):
        # Получаем имя файла без расширения для названия выходного файла
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
        # Получаем список файлов с нужными расширениями
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
