from flask import Flask, request, jsonify
import os
import uuid
import threading

from consensus import build_consensus 

app = Flask(__name__)


tasks_status = {}

def process_task(task_id, dir_in, dir_out):
    try:
        tasks_status[task_id] = {'status': 'started', 'preparedness': 0, 'message': ''}

        host_fs_prefix = '/data'
        dir_in_container = os.path.join(host_fs_prefix, dir_in.lstrip('/'))
        dir_out_container = os.path.join(host_fs_prefix, dir_out.lstrip('/'))
        print(dir_out_container)

        
        if not os.path.isdir(dir_in_container):
            tasks_status[task_id]['status'] = 'error'
            tasks_status[task_id]['message'] = f'Входная директория {dir_in} не найдена.'
            return

        
        os.makedirs(dir_out_container, exist_ok=True)

        
        adapter_path = os.path.join(dir_in_container, 'adapter.fasta')
        if not os.path.isfile(adapter_path):
            tasks_status[task_id]['status'] = 'error'
            tasks_status[task_id]['message'] = f'Файл адаптера {adapter_path} не найден.'
            return

        
        supported_extensions = ('.fasta', '.fastq', '.sam')
        files = [f for f in os.listdir(dir_in_container) if f.endswith(supported_extensions)]
        if not files:
            tasks_status[task_id]['status'] = 'error'
            tasks_status[task_id]['message'] = f'В директории {dir_in} нет файлов с прочтениями.'
            return

        total_files = len(files)
        processed_files = 0

        dir_out_container += '/'
        
        count_right_files = 0
        for filename in files:
            file_path = os.path.join(dir_in_container, filename)
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            name_consensus = f'{base_name}_consensus.fastq'
            if base_name != 'adapter':
                tmp = build_consensus(
                    adapter_fasta=adapter_path,
                    path_to_reads=file_path,
                    path_to_outdir=dir_out_container,
                    muscle_bin_full_path='/app/muscle3.8.31_i86linux64',
                    name_consensus=name_consensus,
                    path_to_minimap2='/app/minimap2-2.28_x64-linux/minimap2'
                )
                if tmp == 1:
                    count_right_files += 1
            else:
                continue

            processed_files += 1
            
            tasks_status[task_id]['preparedness'] = int((processed_files / (total_files - 1)) * 100)

        
        tasks_status[task_id]['status'] = 'completed'
        tasks_status[task_id]['message'] = f'Обработано {count_right_files} из {total_files - 1} файлов'

    except Exception as e:
        tasks_status[task_id]['status'] = 'error'
        tasks_status[task_id]['message'] = str(e)
        print(e)

@app.route('/api/process', methods=['POST'])
def process_reads():
    data = request.get_json()
    if not data:
        return jsonify({'error': 'Пустое тело запроса или некорректный формат JSON'}), 400

    dir_in = data.get('dir_in')
    dir_out = data.get('dir_out')


    
    task_id = str(uuid.uuid4())

    
    thread = threading.Thread(target=process_task, args=(task_id, dir_in, dir_out))
    thread.start()

    return jsonify({'id': task_id}), 200

@app.route('/api/status', methods=['GET'])
def get_status():
    task_id = request.args.get('id')
    if not task_id:
        return jsonify({'error': 'Необходимо указать параметр id'}), 400

    if task_id not in tasks_status:
        return jsonify({'error': 'Задача с указанным id не найдена'}), 404

    status = tasks_status[task_id]
    return jsonify({
        'id': task_id,
        'status': status.get('status', ''),
        'preparedness': status.get('preparedness', 0),
        'message': status.get('message', '')
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
