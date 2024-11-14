**Инструкция по запуску Docker-образа с REST API для создания консенсусных последовательностей**

Чтобы запустить Docker-образ и обработать полимеразные прочтения для получения консенсусных последовательностей, используйте следующую команду:

```bash
docker run -it --rm -v /путь/к/папке/с/прочтениями:/data -p 5000:5000 divolin/cons_for_oms_app:flask_latest 

```

После этого запустится flask и можно будет делать api запросы локально по запросу http://localhost:5000





**API запросы:**


**Отправка на обработку папку с полимеразными прочтениями:**



**Отправка запроса с `requests`:**

```python
import requests

url = 'http://localhost:5000/api/process_json'
headers = {'Content-Type': 'application/json'}
json = {
    'dir_in': '/path/to/input',
    'dir_out': '/path/to/output'
}

response = requests.post(url, json=json, headers=headers)
print(response.json())
```

**Отправка запроса с cURL:**

```bash
curl -X POST -H "Content-Type: application/json" -d '{"dir_in": "/path/to/input", "dir_out": "/path/to/output"}' "http://localhost:5000/api/process_json"
```

**Получение статуса:**



**Отправка запроса с `requests`:**

```python
import requests

url = 'http://localhost:5000/api/get_status'
params = {'id': '12345'}

response = requests.get(url, params=params)
print(response.json())
```

**Отправка запроса с cURL:**

```bash
curl -X GET "http://localhost:5000/api/get_status?id=12345"
```







