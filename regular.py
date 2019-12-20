#import re
#
#str = "genomic RNAS"
#print(re.findall(r'[genomic]{7}.{1}[RNA]{3}', str))

import json

json_dict = {}

json_dict["NC_000128"] = {}
json_dict["NC_000128"]["name"] = "genomic RNA"
json_dict["NC_000128"]["test"] = {"genomic sNA": 213}

json_dict["NC_000129"] = {"name": "genomic RNA"}
#json_dict["NC_000129"]["name"] = "genomic RNA"
for x in json_dict:
    print(x)

print(json_dict)

with open('data.json', 'w') as f:
    json.dump(json_dict, f)