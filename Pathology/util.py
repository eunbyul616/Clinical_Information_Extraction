import pandas as pd
import numpy as np
import os
import glob
import re
from tqdm import tqdm
import collections
import json
tqdm.pandas()

"""

DICTIONARY_PATH: 사전 파일 위치
organ_fname: organ명 사전 경로
opname_fnmae: 시술명 사전 경로
lymph_node_fname: 림프절 코드 전처리 사전 (e.g. peribronchial ln -> #13)

"""

DICTIONARY_PATH = './dict/'

organ_fname = 'organ.txt'
opname_fname = 'opname.txt'
lymph_node_fname = 'lymph_node_dict.json'

with open(os.path.join(DICTIONARY_PATH, organ_fname)) as f:
    organ = [line.strip() for line in f.readlines() if len(line.strip()) > 0]

with open(os.path.join(DICTIONARY_PATH, opname_fname)) as f:
    opname = [line.strip() for line in f.readlines() if len(line.strip()) > 0]

organ_pattern = '|'.join(['('+i+')' for i in organ])
opname_pattern = '|'.join(['('+i+')' for i in opname])

with open(os.path.join(DICTIONARY_PATH, lymph_node_fname)) as f:
    ln_dict = json.load(f)

def extract_diagnosis(x):
    """

    extract diagnosis from full description

    """

    result = []
    main_script = x
    
    x = re.split('diagnosis[\s]*:', x.lower())[1:]
    x = [d for d in x if len(d.strip()) > 0]
   
    if len(x) == 0:
        x = main_script
        x = [x.lower()]
    
    for diagnosis in x:
        try:
            diagnosis = re.split('(note[\s]*:)|(gross[\s:]*)', diagnosis)[0]
            diagnosis = re.sub('(?P<bracket1>\()(?P<char>.+?)(?P<bracket2>\))', '( \g<char> )', diagnosis)
            diagnosis = re.sub('[\s]{2,}', ' ', diagnosis).strip()
            diagnosis = re.sub('[-]{2,}', '', diagnosis).strip()

            list_ = []
            for idx, w in enumerate(diagnosis.strip().split()):
                if idx < len(diagnosis.strip().split()) - 1:
                    if re.sub(',', '', diagnosis.strip().split()[idx+1].lower()).strip() in organ:
                        list_.append(re.sub('(?P<char>([0-9]+[&,-]*)+\))', 'a)', w))
                    else:
                        list_.append(w)
                else:
                    list_.append(w)

            diagnosis = ' '.join(list_)
            
            # delete Korean
            diagnosis = re.sub('[가-힣]', '', diagnosis)
            
            diagnosis = re.sub('(?P<char>([A-Za-z0-9]+[&,-]*)+\))', '\n\g<char>', diagnosis).strip()
            diagnosis = re.sub('(?P<char>([A-Za-z]+[&,-]*)+\))', '\n\g<char>', diagnosis).strip()
            diagnosis = re.sub('\( see note \)', '', diagnosis).strip()
            
            if re.search('[a-z]\)', diagnosis.strip().split()[0]):
                diagnosis = ' '.join(diagnosis.split()[1:])
                diagnosis = [sent.strip() for sent in re.split('[A-Za-z]\)', diagnosis) if len(sent.strip()) > 0]
                
                result.extend(diagnosis)
                
            else:
                result.append(diagnosis)
                
        except Exception:
            pass
    
    return result

def split_sent(x):
    """

    split first line and histologic diagnosis

    """

    result = []
    first_sent = ''
    histology = ''
    for row in x:
        row = [r.strip() for r in re.split(':', row) if len(r.strip()) > 0]
        if len(row) > 1: 
            first_sent = row[0]
            first_sent = first_sent.replace('(', ',(').replace(')', '),')
            first_sent = re.sub('[\(\)]', '', first_sent)
            first_sent = re.sub('and', ',', first_sent)
            first_sent = [w.strip() for w in re.split('[,]', first_sent) if len(w.strip()) > 0]
            
            try:
                histology = ' '.join(row[1:])
            except Exception:
                histology = ''
                                 
        result.append((first_sent, histology))
        
    return result

def labeling(x):
    """

    extract organ, location, operation name from first line

    """
    result_dict = collections.defaultdict(list)
    
    for i in x:
        i = i.strip()
        i = i.replace('-', ' ')
        i = re.sub('[\s]{2,}', ' ', i)
        
        if re.search(organ_pattern, i):
            result_dict['organ'].append(re.search(organ_pattern, i).group())
        
        location_pattern = '(right)|(left)|(rul)|(rll)|(lul)|(lll)|(#[0-9]+[rl]*)|([#]*[0-9]+[rl]{1})|' \
                            + '|'.join(['({})'.format(i) for i in list(ln_dict.keys())])
        
        if re.search(location_pattern, i):
            if not re.search('[0-9]+( )*th', i):
                result_dict['location'].append(i)
        
        if re.search(opname_pattern, i):
            result_dict['opname'].append(re.search(opname_pattern, i).group())
        
    return result_dict

def covert_list_to_str(x, column_name):
    """

    convert list to str
    e.g. ['left upper lobe', 'left lower lobe'] -> 'left upper lobe, left lower lobe'

    """

    try:
        return ', '.join(set(x[column_name]))
    
    except KeyError:
        return ''
    
def find_lymph_node(x):
    """

    find information about lymph nodes
    e.g. ( right ln 13, 0/2; peribronchial ln, 0/2; ln #2r, 0/1; ln #3a, 0/1; ln #4r, 0/2; ln #7, 0/2; ln #9, 0/1; ln #10, 0/1; ln #11s, 0/1; ln #13s, 0/1 )

    """

    result = []
    x = re.sub('\(\s*[a-z]{1}\s*\)', '', x)
    
    for s in re.findall('\(.+?\)', x):
        if re.search('ln #[0-9]+', s) and (re.search('[0-9]+/[0-9]+', s)):
            result.append(s)
            
    return result

def find_tumor_size(x):
    """

    extract tumor size
    
    e.g. 
        2.4 x 2 x 2 cm -> 2.4
        0.9 cm in greatest dimension -> 0.9

    """

    result = []
    
    try:
        for s in re.findall('[0-9.\s]+x[0-9.\s]+x[0-9.\s]+cm', x):
            s = re.sub('[\s]+', '', s)
            result.append(max([float(i.replace('cm', '').strip()) for i in s.split('x')]))

        for s in re.findall('[0-9.\s]+cm in greatest dimension', x):
            s = re.sub('[\s]+', '', s)
            result.append(re.search('[0-9.\s]+cm', s).group().replace('cm', '').strip())

    except Exception:
        pass
    
    return ', '.join(map(str, result))

def find_direction(x):
    """

    extract direction from 'location' column

    """

    result = []
    
    try:
        x = re.sub('#[\s]*(?P<num>[0-9]+)r', '\g<num>right', x)
        x = re.sub('#[\s]*(?P<num>[0-9]+)l', '\g<num>left', x)
        
        for s in re.findall('(right)|(left)', x):
            result.extend(s)
        result = [i for i in result if len(i) > 1]
        
    except TypeError:
        pass
    
    return list(set(result))

def preprocess_ln(x):
    """

    preprocessing informaiton about lymph nodes
    e.g. #12u -> #12
        
    """
    result = []
    try:
        for j in x:
            for k, v in ln_dict.items():
                j = re.sub(k, v, j)
            result.extend([i.replace('ln', '').strip() for i in re.split('[;.,]', re.sub('[\(\)]','' , j))])
            
    except Exception:
        pass
    
    return result