from flask import Flask, request, jsonify
import os
import uuid
import threading
import subprocess
import traceback
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from consensus import build_consensus, make_accuracy  

app = Flask(__name__)

tasks_status = {}
accuracy_tasks_status = {}


def process_file(filename, dir_in_container, dir_out_container, adapter_path):
    base_name = os.path.splitext(os.path.basename(filename))[0]
    name_consensus = f'{base_name}_consensus.fastq'
    if base_name == 'adapter':
        return 0  
    tmp = build_consensus(
        adapter_fasta=adapter_path,
        path_to_reads=filename,
        path_to_outdir=dir_out_container,
        muscle_bin_full_path='/app/muscle3.8.31_i86linux64',
        name_consensus=name_consensus,
        path_to_minimap2='/app/minimap2-2.28_x64-linux/minimap2'
    )
    return 1 if tmp == 1 else 0


def process_fasta_file(file_path, output, folder_name, subread_prefix):
    """
    Обрабатывает файл fasta и записывает прочтения в выходной файл.
    """
    with open(file_path, 'r') as input_file:
        read_number = 1
        sequence_name = None
        sequence_data = []
        for line in input_file:
            line = line.strip()
            if line.startswith('>'):  
                if sequence_name and sequence_data: 
                    output.write(f">{folder_name}:{subread_prefix}{read_number}\n")
                    output.write(''.join(sequence_data) + '\n')
                    read_number += 1
                sequence_name = line
                sequence_data = []
            else:
                sequence_data.append(line) 
        if sequence_name and sequence_data:
            output.write(f">{folder_name}:{subread_prefix}{read_number}\n")
            output.write(''.join(sequence_data) + '\n')


def combine_subreads(input_folder):
    output_file = os.path.join(input_folder, 'subreads.fasta')
    with open(output_file, 'w') as output:
        for root, _, files in sorted(os.walk(input_folder)):
            folder_name = os.path.basename(root)  
            if 'real_reads.fasta' in files:
                process_fasta_file(
                    os.path.join(root, 'real_reads.fasta'), 
                    output, 
                    folder_name, 
                    subread_prefix="subread_"
                )
            if 'real_compl_reads.fasta' in files:
                process_fasta_file(
                    os.path.join(root, 'real_compl_reads.fasta'), 
                    output, 
                    folder_name, 
                    subread_prefix="subread_compl_"
                )


def combine_fastq_files(input_folder, output_file):
    with open(output_file, 'w') as output:
        for filename in sorted(os.listdir(input_folder)):
            if filename.endswith('.fastq'):
                if filename == 'consensuses.fastq':
                    continue
                file_path = os.path.join(input_folder, filename)
                with open(file_path, 'r') as input_file:
                    lines = input_file.readlines()
                    if len(lines[3]) != len(lines[1]):
                        continue 
                    sequence_name = filename.rsplit('.', 1)[0]
                    output.write(f"@{sequence_name}\n")
                    output.write(lines[1])  
                    output.write("+\n")
                    output.write(lines[3]) 


def process_task(task_id, dir_in, dir_out):
    try:
        start_time = time.time()
        # Инициализируем статус задачи
        tasks_status[task_id] = {'status': 'started', 'preparedness': 0, 'message': ''}

        host_fs_prefix = '/data'
        dir_in_container = os.path.join(host_fs_prefix, dir_in.lstrip('/'))
        dir_out_container = os.path.join(host_fs_prefix, dir_out.lstrip('/'))
        print(dir_out_container)

        # Проверяем входную директорию
        if not os.path.isdir(dir_in_container):
            tasks_status[task_id]['status'] = 'error'
            tasks_status[task_id]['message'] = f'Input directory {dir_in} not found.'
            return

        os.makedirs(dir_out_container, exist_ok=True)

        # Проверяем наличие адаптерного файла
        adapter_path = os.path.join(dir_in_container, 'adapter.fasta')
        if not os.path.isfile(adapter_path):
            tasks_status[task_id]['status'] = 'error'
            tasks_status[task_id]['message'] = f'Adapter file {adapter_path} not found.'
            return

        # Фильтруем входные файлы
        supported_extensions = ('.fasta', '.fastq', '.sam')
        files = [os.path.join(dir_in_container, f) for f in os.listdir(dir_in_container) if f.endswith(supported_extensions)]
        files = sorted(files)
        if not files:
            tasks_status[task_id]['status'] = 'error'
            tasks_status[task_id]['message'] = f'No read files found in directory {dir_in}.'
            return

        total_files = len(files)
        dir_out_container += '/'

        count_right_files = 0
        processed_files = 0

        base_names = [os.path.splitext(os.path.basename(f))[0] for f in files]
        has_adapter = 'adapter' in base_names
        effective_total = total_files - 1 if has_adapter else total_files

        with ProcessPoolExecutor(max_workers=8) as executor:
            futures = {
                executor.submit(process_file, f, dir_in_container, dir_out_container, adapter_path): f
                for f in files
            }

            for future in as_completed(futures):
                result = future.result()
                if result == 1:
                    count_right_files += 1
                processed_files += 1

                if effective_total > 0:
                    tasks_status[task_id]['preparedness'] = int(((processed_files-1) / effective_total) * 100)
        with open(os.path.join(dir_out_container, "log.txt"), "a") as f:
            print("--------", time.time() - start_time, "-----------", file=f)
        print("--------", time.time() - start_time, "-----------")
        combine_subreads(dir_out_container)
        combine_fastq_files(dir_out_container, os.path.join(dir_out_container, 'consensuses.fastq'))
        tasks_status[task_id]['status'] = 'completed'
        tasks_status[task_id]['message'] = f'Processed {count_right_files} out of {effective_total} files'

    except Exception as e:
        tasks_status[task_id]['status'] = 'error'
        tasks_status[task_id]['message'] = str(e)



def process_accuracy(task_id, reference, reads):
    """
    Выполняет задачу по вычислению точности для заданных ридов и референса в фоновом потоке.
    Обновляет статус задачи в словаре accuracy_tasks_status.
    После вычисления точности удаляет временные файлы.
    """
    accuracy_tasks_status[task_id] = {'status': 'started', 'preparedness': 0, 'message': ''}

    temp_files = [
        'two_reads.fasta',
        'two_reads_compl.fasta',
        'two_reads_allig.fasta',
        'two_reads_compl_allig.fasta'
    ]

    try:
        host_fs_prefix = '/data'
        reference_path = os.path.join(host_fs_prefix, reference.lstrip('/'))
        reads_path = os.path.join(host_fs_prefix, reads.lstrip('/'))

        if not os.path.isfile(reference_path):
            accuracy_tasks_status[task_id]['status'] = 'error'
            accuracy_tasks_status[task_id]['message'] = f'Reference file {reference} not found.'
            return

        if not os.path.isfile(reads_path):
            accuracy_tasks_status[task_id]['status'] = 'error'
            accuracy_tasks_status[task_id]['message'] = f'Reads file {reads} not found.'
            return
        accuracy = make_accuracy(reference_path, reads_path)

        accuracy_tasks_status[task_id]['status'] = 'completed'
        accuracy_tasks_status[task_id]['preparedness'] = 100
        accuracy_tasks_status[task_id]['message'] = f'Calculated accuracy: {accuracy}%'

    except Exception as e:
        accuracy_tasks_status[task_id]['status'] = 'error'
        accuracy_tasks_status[task_id]['message'] = str(e)
    finally:
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)


@app.route('/api/process', methods=['POST'])
def process_reads():
    """
    Эндпоинт для запуска задачи по построению консенсусов.
    Принимает dir_in, dir_out в JSON.
    Возвращает id задачи.
    """
    data = request.get_json()
    if not data:
        return jsonify({'error': 'Empty request body or invalid JSON format'}), 400

    dir_in = data.get('dir_in')
    dir_out = data.get('dir_out')

    task_id = str(uuid.uuid4())

    thread = threading.Thread(target=process_task, args=(task_id, dir_in, dir_out))
    thread.start()

    return jsonify({'id': task_id}), 200

@app.route('/api/status', methods=['GET'])
def get_status():
    """
    Эндпоинт для получения статуса задачи по построению консенсуса.
    Принимает параметр id.
    Возвращает статус, сообщение и прогресс.
    """
    task_id = request.args.get('id')
    if not task_id:
        return jsonify({'error': 'Parameter id is required'}), 400

    if task_id not in tasks_status:
        return jsonify({'error': 'Task with specified id not found'}), 404

    status = tasks_status[task_id]
    return jsonify({
        'id': task_id,
        'status': status.get('status', ''),
        'preparedness': status.get('preparedness', 0),
        'message': status.get('message', '')
    })

@app.route('/api/accuracy', methods=['POST'])
def start_accuracy():
    """
    Эндпоинт для запуска задачи вычисления точности.
    Принимает reference и reads в JSON.
    Возвращает id задачи.
    """
    data = request.get_json()
    if not data:
        return jsonify({'error': 'Empty request body or invalid JSON format'}), 400

    reference = data.get('reference')
    reads = data.get('reads')

    if not reference or not reads:
        return jsonify({'error': 'Missing "reference" or "reads" parameter'}), 400

    task_id = str(uuid.uuid4())

    thread = threading.Thread(target=process_accuracy, args=(task_id, reference, reads))
    thread.start()

    return jsonify({'id': task_id}), 200

@app.route('/api/accuracy_status', methods=['GET'])
def get_accuracy_status():
    """
    Эндпоинт для получения статуса задачи вычисления точности.
    Принимает параметр id.
    Возвращает статус, сообщение и прогресс.
    """
    task_id = request.args.get('id')
    if not task_id:
        return jsonify({'error': 'Parameter id is required'}), 400

    if task_id not in accuracy_tasks_status:
        return jsonify({'error': 'Task with specified id not found'}), 404

    status = accuracy_tasks_status[task_id]
    return jsonify({
        'id': task_id,
        'status': status.get('status', ''),
        'preparedness': status.get('preparedness', 0),
        'message': status.get('message', '')
    })


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
