# -*- coding: utf-8 -*-

import json
import requests
from pathlib import Path
import datetime

# Welcome message
print('\n=== Figshare file downloader V1.0 ===')
print('=== Figshare文件下载工具 V1.0        ===')
print('===  作者：朱云涛      by yuntao Zhu ===\n')


# Set the base URL and ITEM_ID
BASE_URL = 'https://api.figshare.com/v2'
print('在下面输入figshare的文章ID号')
ITEM_ID = input('Enter the item ID from figshare (example:19179137): ')

#Retrieve public metadata from the endpoint
r=requests.get(BASE_URL + '/articles/' + str(ITEM_ID))
#Load the metadata as JSON
if r.status_code != 200:
    print('Something is wrong:',r.content)
else:
    metadata=json.loads(r.text)
    print('Metadata for item ID', ITEM_ID, 'retrieved successfully.')
    # Save metadata to a file
    with open('metadata_' + str(ITEM_ID) + '.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    print('Metadata saved to metadata_' + str(ITEM_ID) + '.json')
    print('元数据已经检索和保存完成！\n')

# Authoriazation token for API calls
api_call_headers = {'Authorization': 'token dkd8rskjdkfiwi49hgkw...'} 
# You may use your own token, just change the value after token
# example: {'Authorization': 'token dkd8rskjdkfiwi49hgkw...'}
# To get your own token ,see https://help.figshare.com/article/how-to-get-a-personal-token
# 请用自己的token，怎么申请参照https://help.figshare.com/article/how-to-get-a-personal-token

# Or use this test set of ids that have small files (To use, delete the '#' in the next line)
item_ids = [ITEM_ID]

file_info = [] #a blank list to hold all the file metadata

print('Downloading files for item IDs:', item_ids)

for i in item_ids:
    r = requests.get(BASE_URL + '/articles/' + str(i) + '/files')
    file_metadata = json.loads(r.text)
    for j in file_metadata: #add the item id to each file record- this is used later to name a folder to save the file to
        j['item_id'] = i
        file_info.append(j) #Add the file metadata to the list

#Download each file to a subfolder named for the article id and save with the file name
for k in file_info:
    response = requests.get(BASE_URL + '/file/download/' + str(k['id']), headers=api_call_headers)
    Path(str(k['item_id'])).mkdir(exist_ok=True)
    open(str(k['item_id']) + '/' + k['name'], 'wb').write(response.content)
    
print('All done! Files downloaded to folders named by item ID.')
print('所有文件已下载到以文章ID命名的文件夹中')
print('=== Figshare file downloader V1.0 ===\n')
